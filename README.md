# ifcbproc
Tool for processing data from the IFCB

## Example usage
The CLI is designed to be extremely easy to script for. If you have the typical IFCB file structure, for example:
```
JC278/
    D20250529/
        D20250529T092035_IFCB225.adc
        D20250529T092035_IFCB225.hdr
        D20250529T092035_IFCB225.roi
        D20250529T094145_IFCB225.adc
        D20250529T094145_IFCB225.hdr
        (...)
    D20250530/
        (...)
    (...)
```
You could use the following bash script:
```bash
for folder in JC278/* ; do
    dayfolder=$(basename $folder)
    echo "Processing $dayfolder"
    python3 ~/ifcbproc/cli.py ecotaxa ./JC278/$dayfolder/*.roi -o ./JC278_EcoTaxa/$dayfolder.zip --maxsize 64m &
done
wait
```
to create an EcoTaxa compatible zip file for every date, split into zip files of around 64MB each for easy uploading to EcoTaxa:
```
JC278_EcoTaxa/
    D20250529_part1.zip
    D20250529_part2.zip
    D20250530_part1.zip
    D20250530_part2.zip
    D20250530_part3.zip
    D20250530_part4.zip
    D20250531_part1.zip
    D20250531_part2.zip
    (...)
```
