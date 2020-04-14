'''
从 mosdepth 结果里整理每个碱基深度，然后就可以用 R 脚本画深度图。
默认输出文件格式为 csv 。
'''

import pathlib
import sys
import gzip

print("python ExtractDepth.py Input.per-base.bed.gz Output.csv\n")

argvs = sys.argv
in_path = pathlib.Path(argvs[1])
out_path = pathlib.Path(argvs[2])
print(in_path)
print(out_path)

with gzip.open(in_path, mode="rb") as i, open(out_path, 'w') as o:
    o.write('Position\tDepth\n')
    for line in i:
        if line.strip():
            line2 = line.decode(encoding='utf-8').strip()
            items = line2.strip().split('\t')
            print(items)
            p1 = int(items[1])
            p2 = int(items[2]) - 1
            d = items[3]
            o.write(f'{p1},{d},{p2},{d}\n')

print('完成')