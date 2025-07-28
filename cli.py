import numpy
import pickle
import pyarrow.parquet as paparquet
import pyarrow
from typing_extensions import Annotated
from libifcb import ROIReader
import time
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
    dt = datetime.strptime(id_split[0], "D%Y%m%dT%H%M%S")
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

def to_ecotaxa(roi_bin_list, out_file, verbose = False, no_image = False, max_size = None):

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

    samples = []
    for roi_bin in roi_bin_list:
        sample = ROIReader(roi_bin[0], roi_bin[1], roi_bin[2])
        if len(sample.rois) > 0:
            samples.append((sample, roi_bin))

            for key in sample.rois[0].trigger.raw:
                valtype = "f"
                if re.match(r'^-?\d+(?:\.\d+)$', sample.rois[0].trigger.raw[key]) is None:
                    valtype = "t"
                adc_keys.add(("acq_" + key, valtype))


    for adc_key in adc_keys:
        ecotaxa_mapping_order.append(adc_key[0])
        ecotaxa_type_def.append(adc_key[1])

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
    join_srcs = []
    join_dsts = []
    join_files = []
    autoname = False
    verbose = False
    recurse = False
    file_heap = []
    for arg in eargs:
        if arg.startswith("--"):
            if arg == "--help":
                command = "help"
                ehelp_msg = None
                help_flag = True
            elif arg == "--output":
                mode_stack.append(mode)
                mode = "output_capture"
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
            elif arg == "--autoname":
                autoname = True
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
            elif mode == "join_capture_src":
                join_srcs.append(arg)
                mode = "join_capture_dst"
            elif mode == "join_capture_dst":
                join_dsts.append(arg)
                mode = "join_capture_fn"
            elif mode == "join_capture_fn":
                join_files.append(arg)
                mode = mode_stack.pop()

    if command is None:
        help_flag = True

    if output_file is None:
        if command == "parquet" or command == "ecotaxa" or (command == "features" and autoname == False):
            ehelp_msg = "Missing output path"
            help_flag = True

    if help_flag:
        print("")
        print("ifcbproc - A tool for processing IFCB data")
        print("")
        if ehelp_msg is not None:
            print("ERROR")
            print(ehelp_msg)
            print("")
        print("Common usage:")
        print("    ifcbproc parquet <roi_file> [roi_file...] -o <output_path>")
        print("    ifcbproc ecotaxa <roi_file> [roi_file...] -o <output_zip_file>")
        print("    ifcbproc features --autoname <roi_file> [roi_file...]")
        print("")
    else:
        if command == "parquet" or command == "ecotaxa" or command == "features":

            roi_bin_list = []

            cleaved_file_heap = set()

            for filen in file_heap:
                if filen.endswith(".roi") or filen.endswith(".adc") or filen.endswith(".hdr"):
                    cleaved_file_heap.add(filen[:-4])
                else:
                    cleaved_file_heap.add(filen)

            for filen in cleaved_file_heap:
                roi_bin_list.append([filen + ".hdr", filen + ".adc", filen + ".roi"])

            if command == "parquet":
                to_parquet(roi_bin_list, output_file)
            elif command == "ecotaxa":
                to_ecotaxa(roi_bin_list, output_file, verbose = verbose, no_image = no_image, max_size = max_size)
