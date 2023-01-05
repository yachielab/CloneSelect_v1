# Image processing utils
Python script `process_images.py` provides a batch processing of microscope images. 

## Usage `convert_tif_to_jpg` command
Convert tiff format images to jpg format using ImageMagick `magick` and `convert` command line tools. The ImageMagick should be installed. 

```sh
python process_images.py convert_tif_to_jpg \ 
--image-dir=Path to image directory \ 
--outdir=Path to output image directory 
```

## Usage `merge_two_JPGs` command
Merge two corresponding jpg images from the different channels to generate a composite jpg image.

```sh
python process_images.py merge_two_jpgs \ 
--image-dir=Path to image directory \ 
--outdir=Path to output image directory 
```

