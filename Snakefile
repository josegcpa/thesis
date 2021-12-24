from glob import glob
import os

all_folders = [x.strip(os.sep).split(os.sep)[-1] for x in glob('gfx/*/')]
all_output_pngs = []
for folder in all_folders:
    for F in glob('gfx/{}/*pdf'.format(folder)):
        root = F.split(os.sep)[-1][:-4]
        all_output_pngs.append("gfx_png/{}/{}.png".format(folder,root))

try: os.makedirs('gfx_png')
except: pass
for d in all_folders:
    try: os.makedirs('gfx_png/{}'.format(d))
    except: pass

rule all:
    input:
        all_output_pngs

rule convert_to_png:
    input:
        "gfx/{folder}/{root}.pdf"
    output:
        "gfx_png/{folder}/{root}.png"
    params:
        dpi = "300"
    shell:
        """
        convert -verbose -density {params.dpi} -trim "{input}" -quality 100 -flatten -colorspace sRGB "{output}"
        """