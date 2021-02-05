# @Date: 2020/6/22
# @Author: Matthew

'''
从 Bismark 的 M-bias 结果里整理数据用于 M-bias 画图。
默认输出文件格式为 csv 。
'''

import pathlib
import re
import argparse

parser = argparse.ArgumentParser(description="从 Bismark 的 M-bias 结果里整理数据用于 M-bias 画图。")
parser.add_argument("-i", "--input", dest="IN", help="输入")
parser.add_argument("-o", "--output", dest="OUT", help="输出路径，为 csv 格式")
argvs = parser.parse_args()

in_path = pathlib.Path(argvs.IN).resolve()
out_path = pathlib.Path(argvs.OUT).resolve()
print("输入：{}".format(in_path))
print("输出：{}".format(out_path))

first_line = re.compile("(C\w{2}) context \(R(\d)\)")
with open(in_path, 'r') as i, open(out_path, 'w') as o:
    o.write("Read,Position,Context,Total_Count,Methylated_Count,Unmethylated_Count,Methylated_Percent\n")
    for line in i:
        re_result = first_line.match(line)
        if re_result:
            context = re_result.group(1)
            read_num = re_result.group(2)
        elif line.startswith("==") or line.startswith("position"):
            pass
        elif re.match("^\d", line):
            items = line.strip().split("\t")
            o.write(f"Read{read_num},{items[0]},{context},{items[4]},{items[1]},{items[2]},{items[3]}\n")
        else:
            pass

print("\n完成\n")
