# Author: Matthew

'''
从 mosdepth 结果里整理每个碱基深度，然后就可以用 R 脚本画深度图。
默认输出文件格式为 csv 。
'''

import pathlib
import gzip
import argparse

parser = argparse.ArgumentParser(description="从 mosdepth 结果里整理每个碱基深度")
parser.add_argument("-i", "--input", dest="IN", help="输入文件，为 mosdepth 的碱基深度结果 *per-base.bed.gz")
parser.add_argument("-o", "--output", dest="OUT", help="结果输出路径，为 csv 格式")
argvs = parser.parse_args()

in_path = pathlib.Path(argvs.IN).resolve()
out_path = pathlib.Path(argvs.OUT).resolve()
print("输入：{}".format(in_path))
print("输出：{}".format(out_path))

with gzip.open(in_path, mode="rb") as i, open(out_path, 'w') as o:
    o.write('Position,Depth\n')
    for line in i:
        if line.strip():
            line2 = line.decode(encoding='utf-8').strip()
            items = line2.strip().split('\t')
            p1 = int(items[1])
            p2 = int(items[2]) - 1
            d = items[3]
            o.write(f'{p1},{d}\n{p2},{d}\n')

print('\n完成\n')