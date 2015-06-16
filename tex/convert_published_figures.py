import os
import subprocess

with open('published_figures','r') as f:
    figures = [x.strip() for x in f.readlines()]

for row in figures:
    print(row)
    if row[-4:] == '.pdf':
        print("...converting")
        result = subprocess.call(['convert', '-density', '150', '-trim',
                                  '-quality', '100', row,
                                  os.path.splitext(row)[0] + ".png"])
        assert result == 0
