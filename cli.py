#!/bin/python3

# Copyright 2025, A Baldwin <alewin@noc.ac.uk>, National Oceanography Centre.
#
# This file is part of ifcbproc, a tool for processing IFCB data.
#
# ifcbproc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License (version 3 only)
# as published by the Free Software Foundation.
#
# ifcbproc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License v3
# along with ifcbproc.  If not, see <http://www.gnu.org/licenses/>.

import numpy
import pickle
import pyarrow.parquet as paparquet
import pyarrow
from typing_extensions import Annotated
from libifcb import ROIReader
import time
import pytz
import csv
import zipfile
import math
import io
from datetime import datetime
import hashlib
import os
import multiprocessing
from multiprocessing import Process
from multiprocessing import shared_memory
from typing import List
import planktofeatures.extractors
import sys
import re


def print_progress(iteration, total, orig_time):
    bar_w = 32
    bar_x = round(bar_w * (iteration/total))
    bar_l = "#"*bar_x
    bar_r = "_"*(bar_w - bar_x)
    if iteration > 1:
        ttime = ((time.time() - orig_time) / iteration) * total
        secs = round(ttime * (1 - (iteration/total)))
        timestr = f"about {secs}s remaining..."
        if secs < 3:
            timestr = "only a few seconds remaining..."
        if secs > 60:
            mins = math.floor(secs / 60)
            secs = secs - (mins * 60)
            timestr = f"about {mins}min {secs}s remaining..."
        if secs > 3600:
            hrs = math.floor(mins / 60)
            mins = mins - (hrs * 60)
            timestr = f"about {hrs}hr {mins}min remaining..."

        print(f"\r[{bar_l}{bar_r}] {timestr}                ", end="")
    else:
        print(f"\r[{bar_l}{bar_r}] getting ready...", end="")

def ifcb_id_to_udt(ifcb_id):
    id_split = ifcb_id.split("_")
    dt = datetime.strptime(id_split[0], "D%Y%m%dT%H%M%S").replace(tzinfo=pytz.UTC)
    res = int(dt.timestamp())
    udt = "udt__usa_mc_lane_research_laboratories__imaging_flow_cytobot__" + id_split[1].lower() + "__" + str(res)
    if len(id_split) > 2:
        imid = int(id_split[2])
        udt = udt + "__" + str(imid)
    #print(udt)
    #print(small_udt(udt))
    return udt

def small_udt(big_udt):
    if big_udt.startswith("udt__"):
        udt_cmps = big_udt[5:].split("__")
        vendor = udt_cmps[0]
        device = udt_cmps[1]
        vdp = hashlib.sha256((vendor + "__" + device).encode("utf-8")).hexdigest()[0:12]
        device_id = udt_cmps[2].lower()
        did = hashlib.sha256((device_id).encode("utf-8")).hexdigest()[0:16]
        timestamp = int(udt_cmps[3])
        ts = '{:012x}'.format(timestamp)
        small_udt = "udt_" + vdp + "_" + did + "_" + ts
        if len(udt_cmps) > 4:
            imid = udt_cmps[4]
            small_udt = small_udt + "_" + hashlib.sha256((imid).encode("utf-8")).hexdigest()[0:16]
        return small_udt
    else:
        return big_udt

def to_parquet(roi_bin_list, out_file):
    schema = pyarrow.schema([
        ("roi_index", pyarrow.int64()),
        ("trigger_number", pyarrow.int64()),
        ("image", pyarrow.struct([
                ("bytes", pyarrow.binary()),
                ("path", pyarrow.string())
            ]))
    ])
    parquet_writer = paparquet.ParquetWriter(out_file, schema)

    for roi_bin in roi_bin_list:
        sample = ROIReader(roi_bin[0], roi_bin[1], roi_bin[2])
        print(len(sample.rois))

        roi_indicies = []
        roi_images = []
        i = 1
        for roi in sample.rois:
            imbuffer = io.BytesIO()
            roi.image.save(imbuffer, "png")
            imbytes = imbuffer.getvalue()
            roi_images.append({"bytes": imbytes, "path": "roi.png"})
            roi_indicies.append(i)
            i += 1



        batch = pyarrow.RecordBatch.from_pydict({
            "roi_index": roi_indicies,
            "trigger_number": roi_indicies,
            "image": roi_images
            }, schema=schema)
        parquet_writer.write_batch(batch)

    parquet_writer.close()

def human_to_bytes(human_size):
    human_size.lower()
    byte_size = int(re.sub("[^\\d]", "", human_size))
    if human_size.endswith("k"):
        byte_size = byte_size * 1024
    elif human_size.endswith("m"):
        byte_size = byte_size * (1024 ** 2)
    elif human_size.endswith("g"):
        byte_size = byte_size * (1024 ** 3)
    elif human_size.endswith("t"):
        byte_size = byte_size * (1024 ** 4)

    return byte_size

def to_ecotaxa(roi_bin_list, out_file, verbose = False, no_image = False, max_size = None, table_map = {}, joins = [], feature_files = []):

    ecotaxa_mapping = {
                "img_file_name": {
                        "type": "t"
                    },
                "img_rank": {
                        "type": "f"
                    },
                "object_id": {
                        "type": "t"
                    },
                "acq_id": {
                        "type": "t"
                    },
                "object_udt": {
                        "type": "t"
                    },
                "sample_id": {
                        "type": "t"
                    },
                "sample_udt": {
                        "type": "t"
                    },
            }

    ecotaxa_mapping_order = list(ecotaxa_mapping.keys())


    ecotaxa_type_def = []
    for idx in ecotaxa_mapping_order:
        ecotaxa_type_def.append(ecotaxa_mapping[idx]["type"])

    adc_keys = set()
    feature_keys = set()

    #csv.field_size_limit(sys.maxsize) # Needed for some features!

    samples = []
    for roi_bin in roi_bin_list:
        sample = ROIReader(roi_bin[0], roi_bin[1], roi_bin[2])
        matched_feature_files = []
        feature_map = {}
        sample_bn = os.path.splitext(os.path.basename(roi_bin[0]))[0]

        for fncand in feature_files:
            if sample_bn in fncand: # Detect if matching name
                matched_feature_files.append(fncand)

        for matched_feature_file in matched_feature_files:
            with open(matched_feature_file) as csvfile:
                try:
                    csvreader = csv.DictReader(csvfile)
                    fst = True
                    for row in csvreader:
                        feature_map[row["roi_number"]] = row
                        if fst:
                            for key in row.keys():
                                valtype = "f"
                                if re.match(r'^-?\d+(?:\.\d+)$', row[key]) is None:
                                    valtype = "t"
                                feature_keys.add((key, valtype))
                            fst = False
                except Exception:
                    print("Error reading feature file " + matched_feature_file)


        if len(sample.rois) > 0:
            samples.append((sample, roi_bin, feature_map))

            for key in sample.rois[0].trigger.raw:
                valtype = "f"
                if re.match(r'^-?\d+(?:\.\d+)$', sample.rois[0].trigger.raw[key]) is None:
                    valtype = "t"
                adc_keys.add(("acq_" + key, valtype))


    for adc_key in adc_keys:
        ecotaxa_mapping_order.append(adc_key[0])
        ecotaxa_type_def.append(adc_key[1])

    process_features = ["feature_extractor"]

    for feature_key in feature_keys:
        output_key_name = "object_" + feature_key[0]
        if feature_key[0] in process_features:
            output_key_name = "process_" + feature_key[0]
        ecotaxa_mapping_order.append(output_key_name)
        ecotaxa_type_def.append(feature_key[1])

    etfn = []
    for etf in ecotaxa_type_def:
        etfn.append("[" + etf + "]")
    ecotaxa_type_def = etfn

    zip_files = 1
    out_zip = None

    out_file_se = os.path.splitext(out_file)
    if max_size is not None:
        out_zip = zipfile.ZipFile(out_file_se[0] + "_part" + str(zip_files) + out_file_se[1], 'w')
    else:
        out_zip = zipfile.ZipFile(out_file, 'w')

    container = os.path.basename(out_file_se[0]) + "_part" + str(zip_files)

    ecotaxa_md = io.StringIO()
    ecotaxa_md_writer = csv.writer(ecotaxa_md, quoting=csv.QUOTE_NONNUMERIC, delimiter='\t', lineterminator='\n')
    ecotaxa_md_writer.writerow(ecotaxa_mapping_order)
    ecotaxa_md_writer.writerow(ecotaxa_type_def)

    sample = None

    running_compressed_size = len(str(ecotaxa_mapping_order))
    running_compressed_size += len(str(ecotaxa_type_def))

    while len(samples) > 0:
        sample = samples.pop()
        bn = os.path.splitext(os.path.basename(sample[1][0]))[0]
        sample_udt = small_udt(ifcb_id_to_udt(bn))
        for roi in sample[0].rois:
            observation_id = bn + "_" + str(roi.index).zfill(5)

            output_image_path = observation_id + ".png"
            if not no_image:
                imbuffer = io.BytesIO()
                roi.image.save(imbuffer, "png")
                imbytes = imbuffer.getvalue()
                running_compressed_size += len(imbytes)
                out_zip.writestr(container + "/" + output_image_path, imbytes)

            object_md = {
                    "img_file_name": output_image_path,
                    "img_rank": 0,
                    "acq_id": bn + "_TN" + str(roi.trigger.index),
                    "object_id": observation_id,
                    "object_udt": small_udt(ifcb_id_to_udt(observation_id)),
                    "sample_id": bn,
                    "sample_udt": sample_udt
                }


            for feature_key in feature_keys:
                output_key_name = "object_" + feature_key[0]
                if feature_key[0] in process_features:
                    output_key_name = "process_" + feature_key[0]
                object_md[output_key_name] = ""
                #print(sample[2][str(roi.index)])
                if str(roi.index) in sample[2].keys():
                    if feature_key[0] in sample[2][str(roi.index)].keys():
                        object_md[output_key_name] = sample[2][str(roi.index)][feature_key[0]]
                else:
                    print("MISSING FEATURE DATA FOR ROI " + str(roi.index))


            for key in roi.trigger.raw.keys():
                object_md["acq_" + key] = roi.trigger.raw[key]

            ecotaxa_line = []
            for idx in ecotaxa_mapping_order:
                if idx in object_md.keys():
                    ecotaxa_line.append(object_md[idx])
                else:
                    ecotaxa_line.append("")
            ecotaxa_md_writer.writerow(ecotaxa_line)
            running_compressed_size += len(str(ecotaxa_line))

            if max_size is not None:
                if running_compressed_size > max_size:
                    zip_files += 1
                    out_zip.writestr(container + "/ecotaxa.tsv", ecotaxa_md.getvalue())
                    out_zip.close()
                    out_file_se = os.path.splitext(out_file)
                    out_zip = zipfile.ZipFile(out_file_se[0] + "_part" + str(zip_files) + out_file_se[1], 'w')
                    ecotaxa_md = io.StringIO()
                    ecotaxa_md_writer = csv.writer(ecotaxa_md, quoting=csv.QUOTE_NONNUMERIC, delimiter='\t', lineterminator='\n')
                    ecotaxa_md_writer.writerow(ecotaxa_mapping_order)
                    ecotaxa_md_writer.writerow(ecotaxa_type_def)
                    running_compressed_size = len(str(ecotaxa_mapping_order))
                    running_compressed_size += len(str(ecotaxa_type_def))

    out_zip.writestr(container + "/ecotaxa.tsv", ecotaxa_md.getvalue())
    out_zip.close()


def patch_files(roi_bin_list, environmental_data_files = [], verbose = False):
    env_data = {}

    for environmental_data_file in environmental_data_files:
        with open(environmental_data_file) as csvfile:
            csvreader = csv.DictReader(csvfile)
            for row in csvreader:
                mk = None
                if "bin" in row.keys():
                    mk = "bin"
                if "sample" in row.keys():
                    mk = "sample"
                if "filename" in row.keys():
                    mk = "filename"
                keyn = os.path.splitext(os.path.basename(row[mk]))[0]
                env_data[keyn] = row


    for roi_bin in roi_bin_list:
        sample_bn = os.path.splitext(os.path.basename(roi_bin[0]))[0]
        if sample_bn in env_data.keys():
            gps_coords = None
            if ("latitude" in env_data[sample_bn].keys()) and ("longitude" in env_data[sample_bn].keys()):
                gps_coords = env_data[sample_bn]["latitude"], env_data[sample_bn]["longitude"]
            elif ("lat" in env_data[sample_bn].keys()) and ("long" in env_data[sample_bn].keys()):
                gps_coords = env_data[sample_bn]["lat"], env_data[sample_bn]["long"]
            elif ("lat" in env_data[sample_bn].keys()) and ("lon" in env_data[sample_bn].keys()):
                gps_coords = env_data[sample_bn]["lat"], env_data[sample_bn]["lon"]
            if gps_coords is not None:
                if verbose:
                    print("Patching " + sample_bn + " with GPS coords [" + gps_coords[0] + "," + gps_coords[1] + "]")
                with open(roi_bin[0], "r+") as f:
                    data = f.read()
                    print(data)
                    f.seek(0)
                    #f.write(output)
                    #f.truncate()
            else:
                if verbose:
                    print("No GPS data for " + sample_bn)

def generate_features_one_file(sample, csv_file):
    first = True
    orig_time = time.time()
    feature_extractor = planktofeatures.extractors.WHOIVersion4()
    with open(csv_file, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=",", quotechar="\"", quoting=csv.QUOTE_MINIMAL)
        adc_row_idx = 0
        total_rows = len(sample.rows)
        for adc_row_obj in sample.rows:
            adc_row_idx += 1
            #adc_data_row = sample.adc_data[adc_row_idx-1]
            img = adc_row_obj.image
            if img is not None:
                feature_object = feature_extractor.process(img)
                row_features = feature_object.values
                if first:
                    first = False
                    akeys = list(map(str.lower,row_features.keys()))
                    csv_writer.writerow(["roi_number", "feature_extractor", *akeys])
                csv_writer.writerow([adc_row_idx, "whoi_v4", *row_features.values()])

            if adc_row_idx % 16 == 0:
                bar_w = 32
                bar_x = round(bar_w * (adc_row_idx/total_rows))
                bar_l = "#"*bar_x
                bar_r = "_"*(bar_w - bar_x)
                ttime = ((time.time() - orig_time) / adc_row_idx) * total_rows
                if adc_row_idx > 5:
                    secs = round(ttime * (1 - (adc_row_idx/total_rows)))
                    timestr = f"about {secs}s remaining..."
                    if secs < 3:
                        timestr = "only a few seconds remaining..."
                    if secs > 60:
                        mins = math.floor(secs / 60)
                        secs = secs - (mins * 60)
                        timestr = f"about {mins}min {secs}s remaining..."
                    if secs > 3600:
                        hrs = math.floor(mins / 60)
                        mins = mins - (hrs * 60)
                        timestr = f"about {hrs}hr {mins}min remaining..."

                    print(f"\r[{bar_l}{bar_r}] {timestr}                ", end="")
                else:
                    print(f"\r[{bar_l}{bar_r}] getting ready...", end="")

def generate_features(roi_bin_list, out_folder = None, verbose = False):
    roipaths = []
    for roi_bin_group in roi_bin_list:
        roipaths.append(os.path.dirname(os.path.abspath(roi_bin_group[0])))
    common_path = os.path.commonpath(roipaths)
    #print("Common path " + str(common_path))
    if out_folder is not None:
        common_path = os.path.abspath(out_folder)
        #print("Output path " + str(out_path))

    for roi_bin_group in roi_bin_list:
        basename = os.path.splitext(os.path.basename(roi_bin_group[0]))[0]
        roi_reader = ROIReader(roi_bin_group[0], roi_bin_group[1], roi_bin_group[2])
        output_csv = os.path.join(common_path,basename + "_features.csv")
        print(roi_bin_group[2] + " => " + output_csv)
        generate_features_one_file(roi_reader, output_csv)

if __name__ == "__main__":
    eargs = []
    help_flag = False
    ehelp_msg = "No command specified"
    for arg in sys.argv[1:]:
        default = True
        if len(arg) > 1:
            if arg[0] == "-" and (not arg[1] == "-"):
                default = False
                for c in arg[1:]:
                    if c == "h":
                        eargs.append("--help")
                    elif c == "o":
                        eargs.append("--output")
                    elif c == "r":
                        eargs.append("--recurse")
                    elif c == "v":
                        eargs.append("--verbose")
                    else:
                        ehelp_msg = "Unrecognised flag \"-" + c + "\""
                        help_flag = True
        if default:
            eargs.append(arg)

    mode = "command"
    mode_stack = []
    command = None
    output_file = None
    no_image = False
    max_size = None
    joins = []
    tables = []
    verbose = False
    recurse = False
    multiple_capture_switch = False
    file_heap = []
    for arg in eargs:
        if arg.startswith("--"):
            if multiple_capture_switch:
                multiple_capture_switch = False
                mode = mode_stack.pop() # Break out of the current multiple capture
            if arg == "--help":
                command = "help"
                ehelp_msg = None
                help_flag = True
            elif arg == "--output":
                mode_stack.append(mode)
                mode = "output_capture"
            elif arg == "--table":
                mode_stack.append(mode)
                mode = "table_capture"
            elif arg == "--tables":
                mode_stack.append(mode)
                mode = "tables_capture"
                multiple_capture_switch = True
            elif arg == "--join":
                mode_stack.append(mode)
                mode = "join_capture"
            elif arg == "--maxsize":
                mode_stack.append(mode)
                mode = "maxsize_capture"
            elif arg == "--verbose":
                verbose = True
            elif arg == "--noimage":
                no_image = True
            elif arg == "--recurse":
                recurse = True
            else:
                ehelp_msg = "Unrecognised option \"" + arg + "\""
                help_flag = True
        else:
            if mode == "command":
                if arg == "help":
                    command = "help"
                    ehelp_msg = None
                    help_flag = True
                elif arg == "parquet":
                    command = "parquet"
                    mode = "files"
                elif arg == "ecotaxa":
                    command = "ecotaxa"
                    mode = "files"
                elif arg == "patch":
                    command = "patch"
                    mode = "files"
                elif arg == "features":
                    command = "features"
                    mode = "files"
                else:
                    ehelp_msg = "Unrecognised command \"" + arg + "\""
                    help_flag = True
            elif mode == "files":
                file_heap.append(arg)
            elif mode == "output_capture":
                output_file = arg
                mode = mode_stack.pop()
            elif mode == "maxsize_capture":
                max_size = human_to_bytes(arg)
                mode = mode_stack.pop()
            elif mode == "join_capture":
                joins.append(arg)
                mode = mode_stack.pop()
            elif mode == "table_capture":
                tables.append(arg)
                mode = mode_stack.pop()
            elif mode == "tables_capture":
                tables.append(arg)

    if command is None:
        help_flag = True

    if output_file is None:
        if command == "parquet" or command == "ecotaxa":
            ehelp_msg = "Missing output path"
            help_flag = True

    if help_flag:
        print("")
        print("ifcbproc - A tool for processing IFCB data")
        print("")
        print("Copyright 2025, A Baldwin <alewin@noc.ac.uk>, National Oceanography Centre")
        print("This program comes with ABSOLUTELY NO WARRANTY. This is free software,")
        print("and you are welcome to redistribute it under the conditions of the")
        print("GPL version 3 license.")
        print("")
        if ehelp_msg is not None:
            print("ERROR")
            print(ehelp_msg)
            print("")
        print("Common usage:")
        print("    ifcbproc parquet <roi_file> [roi_file...] -o <output_path>")
        print("    ifcbproc ecotaxa <roi_file> [roi_file...] -o <output_zip_file>")
        print("    ifcbproc features <roi_file> [roi_file...] [-o <output_path>]")
        print("")
    else:
        if command == "parquet" or command == "ecotaxa" or command == "features" or command == "patch":

            roi_bin_list = []

            cleaved_file_heap = set()
            csv_heap = set()

            for filen in file_heap:
                if filen.endswith(".roi") or filen.endswith(".adc") or filen.endswith(".hdr"):
                    cleaved_file_heap.add(filen[:-4])
                elif filen.endswith(".csv"):
                    csv_heap.add(filen)
                else:
                    cleaved_file_heap.add(filen)

            table_map = {}
            for table_file in tables:
                table_basename = os.path.splitext(os.path.basename(table_file))
                table_map[table_basename] = []
                with open(table_file) as csvfile:
                    csvreader = csv.DictReader(csvfile)
                    for row in csvreader:
                        table_map[table_basename].append(row)


            for filen in cleaved_file_heap:
                roi_bin_list.append([filen + ".hdr", filen + ".adc", filen + ".roi"])

            if command == "parquet":
                to_parquet(roi_bin_list, output_file, verbose = verbose)
            elif command == "ecotaxa":
                to_ecotaxa(roi_bin_list, output_file, verbose = verbose, no_image = no_image, max_size = max_size, table_map = table_map, joins = joins, feature_files = csv_heap)
            elif command == "patch":
                patch_files(roi_bin_list, environmental_data_files = tables, verbose = verbose)
            elif command == "features":
                generate_features(roi_bin_list, output_file, verbose = verbose)
            else:
                print(command + " unimplemented!")
