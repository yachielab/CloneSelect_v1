import sys
import os.path
import subprocess
import glob
import re
from commandr import command, Run

CONVERT = '/opt/homebrew/bin/convert'
MAGICK = '/opt/homebrew/bin/magick'

@command('convert_tif_to_jpg')
def convert_TIF_to_JPG(image_dir, outdir):
    i = 0
    for subdir in sorted(glob.glob(image_dir + "*")):
        tif_CH1_path, tif_CH3_path = '', ''
        for tif_file in sorted(glob.glob(subdir + "/*.tif")):
            if 'CH1' in tif_file: # Need to change 
                tif_CH1_path = tif_file
            elif 'CH4' in tif_file: # CR4
                tif_CH3_path = tif_file
        if tif_CH3_path == '' or tif_CH1_path == '':
            continue
        
        # GFP channel
        out_CH1 = os.path.basename(tif_CH1_path)
        out_CH1 = 'conv_%d_%s' % (i, out_CH1.replace('tif', 'jpg'))
        out_CH1 = os.path.join(outdir, out_CH1)
        
        # Hoechst channel
        out_CH3 = os.path.basename(tif_CH3_path)
        out_CH3 = 'conv_%d_%s' % (i, out_CH3.replace('tif', 'jpg'))
        out_CH3 = os.path.join(outdir, out_CH3)

        # Need to be adjusted for each experiment
        # optimized for CA062 C-to-T sample
        CMD_Green = '%s +delete -level 1%%,70%%,0.8 -brightness-contrast 20x40 -channel Green -quality 100%% -crop 720x720+0+0 \"%s\" \"%s\"' % (CONVERT, tif_CH1_path, out_CH1)
        CMD_Blue = '%s +delete -level 5%%,90%%,0.8 -quality 100%% -crop 720x720+0+0 \"%s\" \"%s\"' % (CONVERT, tif_CH3_path, out_CH3)
        
        i += 1 # file index
        subprocess.call(CMD_Green, shell=True)
        subprocess.call(CMD_Blue, shell=True)

@command('merge_two_jpgs')
def merge_two_JPGs(target_dir, outdir):
    images = glob.glob(target_dir + '/*.jpg')
    queries = []
    for image in images:
        image_filename = image.split('/')[-1]
        query_index = image_filename.split('_')[1] #unique id for CH1 and CH3 channels
        queries.append(query_index)

    pair = {}
    for q in queries:
        paths = []
        for image in images:
            ret = re.search(r'conv_%s_' % (q), image)
            if ret:
                paths.append(image)
        pair[q] = sorted(paths)

    for k, channels in sorted(pair.items()):
        outfile = os.path.join(outdir, 'merged_' + channels[0]) #.split('/')[2])
        outfile = os.path.join(outdir, channels[0].split('/')[-1].replace('CH1', 'merged'))
        CMD = '\"%s\" \"%s\" -quality 100%% \"%s\" -gravity center -compose Plus -composite \"%s\"' % (MAGICK, channels[0], channels[1], outfile)
        subprocess.call(CMD, shell=True)

if __name__ == '__main__':
    Run()

    #convert_TIF_to_JPG("/Volumes/SOH_32GB/Keyence/02062022/tiff_images/")
    #merge_two_JPGs("/Volumes/SOH_32GB/Keyence/02062022/")

