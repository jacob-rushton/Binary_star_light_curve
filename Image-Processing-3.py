#!/usr/bin/env python
# coding: utf-8

# In[31]:


from pathlib import Path
import os, sys

from astropy.nddata import CCDData
from astropy.io import fits
from astropy.stats import mad_std

import ccdproc as ccdp
import matplotlib.pyplot as plt
import numpy as np

from convenience_functions import show_image

import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from photutils import make_source_mask
from astropy.stats import sigma_clipped_stats

from photutils import aperture_photometry
from photutils import CircularAnnulus
from photutils.aperture import ApertureStats
from photutils.aperture import CircularAperture

from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))


# In[171]:


# cloudy_path = Path('dao_data/images/im_div_mf/cloudy/')
# clouds = ccdp.ImageFileCollection(cloudy_path)
# clouds.summary['date-obs'][-1]


# In[172]:


# Defining file paths for the data bullshit ############
dark_path = Path('plaskett_calibration/darks')
flats_path = Path('plaskett_calibration/flats')
bias_path = Path('plaskett_calibration/bias')
image_path = Path('plaskett_images')

# Now taking the data bullshit and creating an image file collection #########
dark_frames = ccdp.ImageFileCollection(dark_path)
flats_frames = ccdp.ImageFileCollection(flats_path)
bias_frames = ccdp.ImageFileCollection(bias_path)
image_frames = ccdp.ImageFileCollection(image_path)


# In[173]:


# Trimming the files
darks = dark_frames.files
flats = flats_frames.files
images = image_frames.files
for i in range(len(images)):
    image = CCDData.read(image_path / images[i], unit='adu')
    image_oc = ccdp.subtract_overscan(image, 
                                         overscan=image[2307:, :], 
                                         overscan_axis=0)
    image_trim = ccdp.trim_image(image_oc[:2307, :])
    image_trim.write(image_path / "images_trimmed/{}".format(images[i]), overwrite=True)
    if i < len(darks):
        dark = CCDData.read(dark_path / darks[i], unit='adu')
        dark_oc = ccdp.subtract_overscan(dark, 
                                         overscan=dark[2307:, :], 
                                         overscan_axis=0)
        dark_trim = ccdp.trim_image(dark_oc[:2307, :])
        dark_trim.write(dark_path / "darks_trimmed/{}".format(darks[i]), overwrite=True)
    if i < len(flats):
        flat = CCDData.read(flats_path / flats[i], unit='adu')
        flat_oc = ccdp.subtract_overscan(flat, 
                                         overscan=flat[2307:, :], 
                                         overscan_axis=0)
        flat_trim = ccdp.trim_image(flat_oc[:2307, :])
        flat_trim.write(flats_path / "flats_trimmed/{}".format(flats[i]), overwrite=True)


# In[174]:


print(images[0])


# In[32]:


# Redefining file paths for the data bullshit ############
dark_path = Path('plaskett_calibration/darks/darks_trimmed')
flats_path = Path('plaskett_calibration/flats/flats_trimmed')
image_path = Path('plaskett_images/images_trimmed')

# Now taking the data bullshit and creating an image file collection #########
dark_frames = ccdp.ImageFileCollection(dark_path)
flats_frames = ccdp.ImageFileCollection(flats_path)
image_frames = ccdp.ImageFileCollection(image_path)


# In[176]:


dark_frames.summary


# In[177]:


flats_frames.summary


# In[178]:


image_frames.summary


# In[179]:


# The Hunt For the Exposure;  it's unimportant at this point

# image_name = 'dao_c182_2022_006213.fits.gz'.rstrip()

# image_path = Path(dark_path) / image_name

# hdu_list = fits.open(image_path)
# hdu_list.info()

# hdu = hdu_list[0]
# hdu.header


# with fits.open(all_flats[9]) as hdu_flats:
#     hdu_flats.verify("fix")
# hdu_flats = fits.open(all_flats[10])
# hdu_flats.info()

# hdu = hdu_flats[0]
# hdu.header['EXPOSURE']
# hdul[0].header


# In[180]:


# METHOD 2: ##########################

# A) Creating Raw Master Dark by combining all dark frames
calibrated_darks = dark_frames.files_filtered(include_path=True)

raw_master_dark = ccdp.combine(calibrated_darks,
                             method='median',
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                             unit='adu'
                            )

# Saving the file (if I run this uncommented after file already created, it gets angry)
raw_master_dark.meta["combined"] = True

raw_master_dark.write('plaskett_calibration/raw_master_dark.fits')


# In[181]:


#print(raw_master_dark)
#print(type(raw_master_dark))
# B1) Subtracting Raw Master Dark from each flat frame:
from astropy.io import fits as pyfits

a_flat = CCDData.read(flats_frames.files_filtered(include_path=True)[0], unit='adu')

show_image(a_flat, cmap='gray')

all_flats = flats_frames.files_filtered(include_path=True)
print(all_flats)
print(type(all_flats))

def subtractImg(fout, fitsToSubtract, filelist): # Slightly modified Boley function
    '''takes(fout,fitsToSubtract,filelist)  '''
    import pdb
    hdulref = fits.open(fitsToSubtract)
    head = hdulref[0].header
    refdata = hdulref[0].data
    hdulref.close()

    #fh = open(filelist, "r") 
    iter = 0
    #for line in fh:
    for line in filelist: # This is basically all I changed
        fpath = line.rstrip()
        hdul = fits.open(fpath)

        head = hdul[0].header
        imgdata = hdul[0].data

        Ndata = imgdata.size

        Nref = refdata.size

        if not Ndata == Nref:
            print("Mismatch between data {} and reference {}".format(Ndata, Nref))
            sys.exit()

        NY, NX = imgdata.shape

        corrected = imgdata * 0
        for j in range(NY):
            for i in range(NX):
                if refdata[j][i] > imgdata[j][i]:
                    corrected[j][i] = 0
                else:
                    corrected[j][i] = imgdata[j][i] - refdata[j][i]

        head.set('comment', 'Subtracted {} from file {}'.format(fitsToSubtract, fpath))
        hdul[0].data = corrected

        hdul.writeto(fout + "-" + repr(iter).zfill(3) + ".fits", overwrite=True)
        hdul.close()
        iter += 1


subtractImg("plaskett_calibration/flats/flats_trimmed/flats_minus_rmd/flats_minus_rmd", 
            'plaskett_calibration/raw_master_dark.fits', 
            all_flats)       


# In[182]:


# B2) Combining flats - rmd into Master Flat:
def inv_median(a):
    return 1 / np.median(a)

flats_rmd_path = flats_path / "flats_minus_rmd"
flats_rmd_frames = ccdp.ImageFileCollection(flats_rmd_path)
flats_minus_rmd = flats_rmd_frames.files_filtered(include_path=True)
master_flat = ccdp.combine(flats_minus_rmd,
                             method='median',scale=inv_median,
                             sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                             sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                             unit='adu'
                            )

# Saving the file (if I run this uncommented after file already created, it gets angry)
master_flat.meta["combined"] = True

master_flat.write('plaskett_calibration/master_flat.fits')


# In[183]:


CCDData.read('plaskett_calibration/master_flat.fits').data.min()


# In[184]:


# C3) Subtracting raw master dark from each science frame:
all_images = image_frames.files_filtered(include_path=True)

subtractImg("plaskett_images/images_trimmed/im_minus_rmd/im_minus_rmd", 
            "plaskett_calibration/raw_master_dark.fits", all_images) 


# In[189]:


# C4) Attempt at dividing im_minus_rmd files by master flat:

def divImg(fout, fitsToDivide, filelist):  # Literally Boley's - func swtiched to /
    '''takes(fout,fitsToSubtract,filelist)  '''

    hdulref = fits.open(fitsToDivide)
    head = hdulref[0].header
    refdata = hdulref[0].data
    hdulref.close()

    #fh = open(filelist, "r")
    iter = 0
    #for line in fh:
    for line in filelist:
        fpath = line.rstrip()
        hdul = fits.open(fpath)

        head = hdul[0].header
        imgdata = hdul[0].data

        Ndata = imgdata.size

        Nref = refdata.size

        if not Ndata == Nref:
            print("Mismatch between data {} and reference {}".format(Ndata, Nref))
            sys.exit()

        NY, NX = imgdata.shape

        corrected = imgdata * 0
        for j in range(NY):
            for i in range(NX):
                if refdata[j][i] == 0: # The most notable change I made
                    corrected[j][i] = 0
                else:
                    corrected[j][i] = imgdata[j][i] / refdata[j][i]

        head.set('comment', 'Divided {} from file {}'.format(fitsToDivide, fpath))
        hdul[0].data = corrected

        hdul.writeto(fout + "-" + repr(iter).zfill(3) + ".fits")
        hdul.close()
        iter += 1
        
im_rmd_path = image_path / "im_minus_rmd"
im_rmd_frames = ccdp.ImageFileCollection(im_rmd_path)
all_im_minus_rmd = im_rmd_frames.files_filtered(include_path=True)

        
divImg("plaskett_images/images_trimmed/processed_images/proc_im",
       "plaskett_calibration/master_flat.fits", all_im_minus_rmd)


# # Data analysis

# ## Background estimation
# 
# All of these notes are copied from https://buildmedia.readthedocs.org/media/pdf/photutils/v0.3/photutils.pdf
# 
# Having an accurate estimate of the background noise is important for determining the significance of source detections and for estimating photometric errors.
# 
# Unfortunately, accurate background and background noise estimation is a difficult task. Further, because astronomical images can cover a wide variety of scenes, there is not a single background estimation method that will always be applicable. Photutils provides tools for estimating the background and background noise in your data, but they will likely require some tweaking to optimize the background estimate for your data.
# 
# If the background level and noise are relatively constant across an image, the simplest way to estimate these values is to derive scalar quantities using simple approximations. Of course, when computing the image statistics one must take into account the astronomical sources present in the images, which add a positive tail to the distribution of pixel intensities. For example, one may consider using the image median as the background level and the image standard deviation as the 1-sigma background noise, but the resulting values are obviously biased by the presence of real sources.
# 
# A slightly better method involves using statistics that are robust against the presence of outliers, such as the biweight location for the background level and biweight midvariance or median absolute deviation (MAD) for the background noise estimation. However, for most astronomical scenes these methods will also be biased by the presence of astro- nomical sources in the image.
# 

# In[38]:


import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
get_ipython().run_line_magic('matplotlib', 'inline')


# In[39]:


rmd = CCDData.read("plaskett_calibration/raw_master_dark.fits")
mf = CCDData.read("plaskett_calibration/master_flat.fits", unit='adu')
exraw = CCDData.read('plaskett_images/dao_c182_2022_004803.fits.gz', unit='adu')
exproc = CCDData.read(image_path / "processed_images/proc_im-000.fits", unit='adu')
fig, ax = plt.subplots(2,2,figsize=(10,20))
show_image(rmd, cmap='Greys_r', percl=99.9, fig=fig, ax=ax[0,0], show_colorbar=False)
ax[0,0].set_title('Raw Master Dark', fontsize=20)
show_image(mf, cmap='Greys_r', percl=99.9, fig=fig, ax=ax[0,1])
ax[0,1].set_title('Master Flat', fontsize=20)
show_image(exraw, cmap='Greys_r', percl=99.9, fig=fig, ax=ax[1,0], show_colorbar=False)
ax[1,0].set_title('Example raw image', fontsize=20)
show_image(exproc, cmap='Greys_r', percl=99.95, fig=fig, ax=ax[1,1], show_colorbar=False)
ax[1,1].set_title('Example processed image', fontsize=20)
norm = ImageNormalize(stretch=SqrtStretch())
for i in range(2):
    for j in range(2):
        ax[i,j].tick_params(labelsize=18)
        ax[i,j].set_xlabel('X pixel', fontsize=20)
        ax[i,j].set_ylabel('Y pixel', fontsize=20)
plt.subplots_adjust(wspace=0.3)
plt.savefig("RMD_MF_Img-Comp.jpg")


# In[6]:


# aimgf2.summary['date-obs'][-1]


# In[40]:


from astropy.stats import sigma_clipped_stats
mean, median, std = sigma_clipped_stats(exproc, sigma=3.0)#, iters=5)
print((mean, median, std))


# In[41]:


from photutils import make_source_mask
from astropy.stats import sigma_clipped_stats
mask = make_source_mask(exproc, nsigma=3, npixels=5, dilate_size=11)
mean, median, std = sigma_clipped_stats(exproc, sigma=4.0, mask=mask)
print((mean, median, std))


# In[42]:


from photutils.detection import DAOStarFinder
pathway = image_path / "processed_images"
imframe = ccdp.ImageFileCollection(pathway)
imframes = imframe.files_filtered(include_path=True)[:10]
counts = []
median_bg = []
std_bg = []
for frame in imframes: 
    mask = make_source_mask(exproc, nsigma=3, npixels=5, dilate_size=11)
    mean, median, std = sigma_clipped_stats(exproc, sigma=4.0, mask=mask)
    median_bg.append(median)
    std_bg.append(std)
    daofind = DAOStarFinder(fwhm=3.0, threshold=500.*std)  
    sources = daofind(exproc - median)  
    for col in sources.colnames:  
        sources[col].info.format = '%.8g'  # for consistent table output
    source = sources[(sources['ycentroid']>=1500)]#&(sources['ycentroid']<=600)]
    counts.append(np.mean(source['peak'].value))
#     print(sources) 


# In[43]:


positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
apertures = CircularAperture(positions, r=50.);
norm = ImageNormalize(stretch=SqrtStretch())
fig, ax = plt.subplots(figsize=(30,15))
#plt.imshow(exproc, cmap='Greys_r', origin='lower', norm=norm,
#           interpolation='nearest')
show_image(exproc, cmap='Greys_r', percl=99.95, fig=fig, ax=ax)
apertures.plot(color='r', lw=1.5, alpha=0.5);
# fig.savefig("Identified_DAO_Stars.jpg")


# In[44]:


xpixels = sources['xcentroid'].value.astype(int)
ypixels = sources['ycentroid'].value.astype(int)
xpixels


# In[45]:


xindex0 = np.digitize(xpixels, bins=np.arange(0,1000, 100))
yindex0 = np.digitize(ypixels, bins=np.arange(0,2000, 100))


# In[ ]:


# Separating Day 1 and Day 2 images
proc_frames = ccdp.ImageFileCollection(image_path / "processed_images", 
                                       keywords='*')

for hdu, fname in proc_frames.hdus(return_fname=True):
    if hdu.header['UT-DATE'] == "2022-02-26":
        hdu.writeto(image_path / "processed_images/day1" / fname)
    else:
         hdu.writeto(image_path / "processed_images/day2" / fname)


# In[46]:


img_path1 = Path(image_path / 'processed_images/day1/')
imgf1 = ccdp.ImageFileCollection(img_path1)
img_path2 = Path(image_path / 'processed_images/day2/')
imgf2 = ccdp.ImageFileCollection(img_path2)
imgf1.summary


# In[31]:


imgf2.summary


# ## Aligned images in AstroImageJ

# In[47]:


aimg_path1 = Path(image_path / 'processed_images/aligned/day1/')
aimgf1 = ccdp.ImageFileCollection(aimg_path1)
aimg_path2 = Path(image_path / 'processed_images/aligned/day2/')
aimgf2 = ccdp.ImageFileCollection(aimg_path2)


# In[48]:


aimgf2.summary


# In[49]:


tobs2 = np.array(aimgf2.summary['date-obs'].value)
hr2 = []
min2 = []
sec2 = []
for t in tobs2:
    hr = t[11:13]
    minute = t[14:16]
    sec = t[17:19]
    hr2.append(hr)
    min2.append(minute)
    sec2.append(sec)
hr2 = np.array(hr2).astype(float)
min2 = np.array(min2).astype(float)
sec2 = np.array(sec2).astype(float)
time2 = (hr2-hr2.min())*60 + min2 + sec2/60


# In[50]:


tobs1 = np.array(aimgf1.summary['date-obs'].value)
hr1 = []
min1 = []
sec1 = []
for t in tobs1:
    hr = t[11:13]
    minute = t[14:16]
    sec = t[17:19]
    hr1.append(hr)
    min1.append(minute)
    sec1.append(sec)
hr1 = np.array(hr1).astype(float)
min1 = np.array(min1).astype(float)
sec1 = np.array(sec1).astype(float)
time1 = (hr1-hr1.min())*60 + min1 + sec1/60


# ## Day 1 analysis:

# In[51]:


aimgf1.summary


# In[52]:


# Initial guess:
target = (614,1187)
comp = [(357,1228), (614,1468), (660,1585)]
ourstars = [(614,1187), (357,1228), (614,1468), (660,1585)]


# In[53]:


exproc = CCDData.read(aimg_path1/aimgf1.files[0], unit='adu')
fig, ax = plt.subplots(figsize=(10,10))
show_image(exproc, cmap='Greys', percl=99.9, fig=fig, ax=ax, show_colorbar=False)
apertures = CircularAperture(target, r=11.);
apertures.plot(color='r', lw=2.5, alpha=1);
apertures = CircularAperture(target, r=12.);
apertures.plot(color='xkcd:black', ls='--', lw=1.5, alpha=1);
apertures = CircularAperture(target, r=19.);
apertures.plot(color='xkcd:black', lw=1.5, ls='--', alpha=1);

rvals = [11, 11, 11, 11]
for r, compstar in enumerate(comp):
    apertures = CircularAperture(compstar, r=rvals[r]);
    apertures.plot(color='xkcd:bright blue', lw=3, alpha=1);
    apertures = CircularAperture(compstar, r=rvals[r]+2);
    apertures.plot(color='xkcd:bright pink', lw=2, ls='--', alpha=1);
    apertures = CircularAperture(compstar, r=rvals[r]+8);
    apertures.plot(color='xkcd:bright pink', lw=2, ls='--', alpha=1);
#plt.xlim(250,750)
#plt.ylim(1100,1700)
ax.tick_params(labelsize=18)
plt.xlabel('X pixel', fontsize=20)
plt.ylabel('Y pixel', fontsize=20)
plt.title("Example Day 1 Aperture Image", fontsize=16)
legend_elements = [Line2D([0], [0], lw=3, color='xkcd:bright blue', ls='-', label='Calib. star\naperture'),
                   Line2D([0], [0], lw=3, color='xkcd:bright pink', ls='--', label='Calib. star\nBG annulus'),
                   Line2D([0], [0], color='r', lw=3, label='Target aperture'),
                  Line2D([0], [0], color='k', lw=3, label='Target BG annulus')]

ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor= (1.02, 0.4), ncol=1,
            borderaxespad=0, frameon=True, handlelength=3, fontsize=16)
plt.savefig("Aperture_Image_Full.jpg")


# In[54]:


# Zoomed in:

fig, ax = plt.subplots(figsize=(10,10))
show_image(exproc, cmap='Greys', percl=99.9, fig=fig, ax=ax, show_colorbar=False)
apertures = CircularAperture(target, r=11.);
apertures.plot(color='r', lw=2.5, alpha=1);
apertures = CircularAperture(target, r=12.);
apertures.plot(color='k', ls='--', lw=2, alpha=1);
apertures = CircularAperture(target, r=19.);
apertures.plot(color='k', lw=2, ls='--', alpha=1);

rvals = [11, 11, 11, 11]
for r, compstar in enumerate(comp):
    apertures = CircularAperture(compstar, r=rvals[r+1]);
    apertures.plot(color='xkcd:bright blue', lw=3, alpha=1);
    apertures = CircularAperture(compstar, r=rvals[r+1]+2);
    apertures.plot(color='xkcd:bright pink', lw=2, ls='--', alpha=1);
    apertures = CircularAperture(compstar, r=rvals[r+1]+8);
    apertures.plot(color='xkcd:bright pink', lw=2, ls='--', alpha=1);
plt.xlim(250,750)
plt.ylim(1095,1705)
ax.tick_params(labelsize=18)
plt.xlabel('X pixel', fontsize=20)
plt.ylabel('Y pixel', fontsize=20)
plt.title("Zoomed in Day 1 Apertures", fontsize=20)
legend_elements = [Line2D([0], [0], lw=3, color='xkcd:bright blue', ls='-', label='Calib. star\naperture'),
                   Line2D([0], [0], lw=3, color='xkcd:bright pink', ls='--', label='Calib. star\nBG annulus'),
                   Line2D([0], [0], color='r', lw=3, label='Target aperture'),
                  Line2D([0], [0], color='k', lw=3, label='Target BG annulus')]

# plt.legend(handles=legend_elements, loc='lower left', bbox_to_anchor= (0.04, 1.01), ncol=2,
#             borderaxespad=0, frameon=True, handlelength=3, fontsize=16)
plt.legend(handles=legend_elements, loc='best', fontsize=16)
plt.savefig("day1_apertures_zoomed.jpg")


# In[55]:


# apertures = CircularAperture(ourstars, r=10)

# # testing this method from a website (don't use this one!)


# annulus_apertures = CircularAnnulus(ourstars, r_in=11., r_out=20.)
# aperstats = ApertureStats(exproc.data, annulus_apertures)
# bkg_mean = aperstats.mean
# print(bkg_mean)  
# phot_table = aperture_photometry(exproc, apertures)
# for col in phot_table.colnames:
#     phot_table[col].info.format = '%.8g'  # for consistent table output
# print(phot_table)
# aperture_area = apertures.area_overlap(exproc)
# total_bkg = bkg_mean * aperture_area
# phot_bkgsub = phot_table['aperture_sum'].value - total_bkg

# phot_table['total_bkg'] = total_bkg
# phot_table['aperture_sum_bkgsub'] = phot_bkgsub
# for col in phot_table.colnames:
#     phot_table[col].info.format = '%.8g'  # for consistent table output
# print(phot_table)


# In[56]:


import pandas as pd
data = pd.DataFrame(columns=['time', 'star', 'aperture_sum', 'total_bkg_local', 
                             'aperture_sum_bkgsub_local', 'aperture_sum_bkgsub_global', 'aperstats_sum'])
apertures = CircularAperture(ourstars, r=11.);
annulus_apertures = CircularAnnulus(ourstars, r_in=12., r_out=19.)
bgmeans = []
bgstds = []

for i in range(len(aimgf1.files)):
    exproc = CCDData.read(aimg_path1/aimgf1.files[i], unit='adu')

    mask = make_source_mask(exproc, nsigma=3, npixels=5, dilate_size=11)
    mean, median, std = sigma_clipped_stats(exproc, sigma=4.0, mask=mask)
    bgmeans.append(mean)
    bgstds.append(std)
    phot_table_local = aperture_photometry(exproc - median, apertures)
    
    aperstats = ApertureStats(exproc.data, annulus_apertures)
    sums = aperstats.sum
    bkg_mean = aperstats.mean
    print(bkg_mean)
    phot_table = aperture_photometry(exproc, apertures)
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output
    aperture_area = apertures.area_overlap(exproc)
    total_bkg = bkg_mean * aperture_area
    phot_bkgsub = phot_table['aperture_sum'].value - total_bkg

    phot_table['total_bkg'] = total_bkg
    phot_table['aperture_sum_bkgsub'] = phot_bkgsub
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output
    
    dat = np.array([time1[i]*np.ones(4), phot_table['id'].value.astype(int), phot_table['aperture_sum'].value.astype(float),
                    phot_table['total_bkg'].value.astype(float), phot_table['aperture_sum_bkgsub'].value.astype(float), 
                    phot_table_local['aperture_sum'].value.astype(float), sums]).T
    data = data.append(pd.DataFrame(dat, columns=data.columns))


# In[57]:


bgmeans
phot_table


# In[58]:


binary = data.loc[data.star==1].copy()
comp1 = data.loc[data.star.isin([2,3,4])].copy()


# In[59]:


plt.scatter(binary.time.values, -2.5*np.log10(binary.aperture_sum_bkgsub_local.astype(float)), s=10)


# In[60]:


import astropy.units as u
Vstar2 = 12.00
uVstar2 = 0.21
pstar2 = 1.5002 * u.mas 
upstar2 = 0.0129 * u.mas
Vstar3 = 14.69
uVstar3 = 0.08
pstar3 = 0.9215 * u.mas
upstar3 = 0.0207 * u.mas
Vstar4 = 13.15
uVstar4 = 0.05
pstar4 = 1.8988
upstar4 = 0.0183

Vstars = np.array([Vstar2, Vstar3, Vstar4])
uVstars = np.array([uVstar2, uVstar3, uVstar4])
#np.concatenate(np.repeat([Vstars], 120, axis=0))


# In[61]:


[(357,1228), (614,1468), (660,1585)]


# In[62]:


comp1['Amag'] = -2.5 * np.log10(comp1.aperture_sum_bkgsub_local.values.astype(float));
comp1['Vmag'] = np.concatenate(np.repeat([Vstars], 120, axis=0));
comp1['uVmag'] = np.concatenate(np.repeat([uVstars], 120, axis=0));
# comp1['Vmag'] = np.concatenate(np.repeat([Vstars], 212, axis=0));
# comp1['uVmag'] = np.concatenate(np.repeat([uVstars], 212, axis=0));


# In[63]:


fig, ax = plt.subplots()
binmags = - 2.5 * np.log10(binary.aperture_sum_bkgsub_local.values.astype(float)/comp1.loc[comp1.star==2].aperture_sum_bkgsub_local.values.astype(float))
mag = binmags + 12.00
plt.scatter(time1, mag)
ax.invert_yaxis()


# In[64]:


comp1["Vmag"]


# In[65]:


test = comp1.iloc[0:3]
np.mean(test.Vmag+test.Amag)
np.std(test.Vmag+test.Amag)


# In[66]:


dtime.iloc[0].Amag


# In[67]:


mjdt1 = np.array(aimgf1.summary['mjd-obs'].value)
dVavgs = []
corrs = []
corr1 = []
corr2 = []
corr3 = []
for i,t in enumerate(time1):
    dtime = comp1.loc[comp1.time.values==t]
    dV = dtime.Vmag.values - (2.5*np.log10(dtime.aperture_sum_bkgsub_local.values.astype(float)))
    dVl = dtime.Vmag.values - dtime.uVmag.values - (2.5*np.log10(dtime.aperture_sum_bkgsub_local.values.astype(float)))
    dVh = dtime.Vmag.values + dtime.uVmag.values - (2.5*np.log10(dtime.aperture_sum_bkgsub_local.values.astype(float)))
    dVavg = np.mean(dV)
    dVavgs.append(dVavg)
    dVhavg = np.mean(dVh) - dVavg
    dVlavg = dVavg - np.mean(dVl)
    binmag = 2.5 * np.log10(binary.loc[binary.time.values==t].aperture_sum_bkgsub_local.values.astype(float))
    #plt.scatter(t/60, binmag+dVavg, color='k', s=10)
    C = np.mean(-1*dtime.Amag + dtime.Vmag)
    corr1.append(dtime.iloc[0].Vmag-dtime.iloc[0].Amag)
    corr2.append(dtime.iloc[1].Vmag-dtime.iloc[1].Amag)
    corr3.append(dtime.iloc[2].Vmag-dtime.iloc[2].Amag)
    dC = np.sqrt(np.sum(dtime.uVmag**2)) / len(dtime.uVmag)
    corrs.append(C)
dVavgs = np.array(dVavgs)
corrs = np.array(corrs)


# In[68]:


corr1 = np.array(corr1)
corr2 = np.array(corr2)
corr3 = np.array(corr3)
dC


# In[69]:


fig, ax = plt.subplots(figsize=(10,9))
binmags = - 2.5 * np.log10(binary.aperture_sum_bkgsub_local.values.astype(float))
maxmag = (binmags+dVavgs).min()
plt.scatter(time1/60, (binmags+corrs), color='k', s=20)
plt.fill_between(time1/60, (binmags+corrs-dC), (binmags+corrs+dC), 
                 zorder=0, color='xkcd:light grey')
plt.xlabel(tobs1[0][:10] + ': ' + 'Hours since '+ tobs1[0][11:] + ' UTC', fontsize=20)
plt.ylabel('Apparent V magnitude', fontsize=20)
ax.tick_params(labelsize=18)
plt.axvline(time1[15]/60, color='r', lw=1, ls='--')
plt.fill_between(np.linspace(time1[13],time1[17],10)/60, 13.9, 15.2, zorder=0, color='r', alpha=0.075)
plt.axvline(149.7/60, color='r', lw=1, ls='--')
plt.axhline((binmags+corrs).min(), zorder=0, color='xkcd:grey', ls='--')
plt.fill_between(np.linspace(time1[60],time1[68],10)/60, 13.9, 15.2, zorder=0, color='r', alpha=0.075)
#plt.axvline(0.16, color='r', lw=1, ls='--')
#plt.fill_between(np.linspace(0.155,0.18,10), 13.9, 15.2, zorder=0, color='r', alpha=0.1)
legend_elements = [Line2D([0], [0], marker='o', color='k', ls='', label='Corrected V\nmagnitude'),
                   Line2D([0], [0], color='xkcd:light grey', lw=10, label='$\sigma(V)$'),
                    Line2D([0], [0], color='r', lw=5, ls='--', label='orbit extremas'),
                  Line2D([0], [0], color='r', lw=10, alpha=0.1, label=r'$\sigma(\rm{extremas})$'),
                  Line2D([0], [0], color='xkcd:grey', lw=5, ls='--', label=r'V$_{\rm{max}}$'+'$\simeq${:.3}'.format((binmags+corrs).min()))]
# ax.legend(handles=legend_elements, loc='best', bbox_to_anchor= (1.01, 0.45), ncol=1,
#             borderaxespad=0, frameon=True, handlelength=3, fontsize=16, markerscale=2)
ax.legend(handles=legend_elements, loc='lower center', ncol=1, frameon=True, 
          handlelength=2, fontsize=16, markerscale=2)
plt.ylim(13.95, 15.15)
ax.invert_yaxis()
ax.set_title('Light curve of eclipsing binary J074658', fontsize=22)
plt.savefig("Binary_Light_Curve.png")
print((binmags + corrs).min())


# In[70]:


(149.7-time1[15])*4/60


# In[71]:


5.6/60


# In[72]:


# fig, ax = plt.subplots(1,3,figsize=(20,6))
fig, ax = plt.subplots(3, 1,figsize=(9, 24))


binmags = - 2.5 * np.log10(binary.aperture_sum_bkgsub_local.values.astype(float))
ax[0].scatter(time1/60, (binmags+corr1), color='k', s=20)
ax[0].fill_between(time1/60, (binmags+corr1-0.21), (binmags+corr1+0.21), 
                 zorder=0, color='xkcd:light grey')
ax[0].axhline((binmags+corr1).min(), zorder=0, color='xkcd:grey', ls='--')

legend_elements = [Line2D([0], [0], marker='o', color='k', ls='', label='Corrected V\nmagnitude'),
                   Line2D([0], [0], color='xkcd:light grey', lw=10, label='$\sigma(V)$'),
                  Line2D([0], [0], color='xkcd:grey', lw=5, ls='--', label=r'V$_{\rm{max}}$'+'$\simeq${:.3}'.format((binmags+corr1).min()))]
ax[0].legend(handles=legend_elements, loc='lower center', ncol=1, frameon=True, 
          handlelength=2, fontsize=16, markerscale=2)

ax[1].scatter(time1/60, (binmags+corr2), color='k', s=20)
ax[1].fill_between(time1/60, (binmags+corr2-0.08), (binmags+corr2+0.08), 
                 zorder=0, color='xkcd:light grey')
ax[1].axhline((binmags+corr2).min(), zorder=0, color='xkcd:grey', ls='--')

legend_elements = [Line2D([0], [0], marker='o', color='k', ls='', label='Corrected V\nmagnitude'),
                   Line2D([0], [0], color='xkcd:light grey', lw=2, label='$\sigma(V)$'),
                  Line2D([0], [0], color='xkcd:grey', lw=2, ls='--', label=r'V$_{\rm{max}}$'+'$\simeq${:.3}'.format((binmags+corr2).min()))]
ax[1].legend(handles=legend_elements, loc='lower center', ncol=1, frameon=True, 
          handlelength=2, fontsize=16, markerscale=2)

ax[2].scatter(time1/60, (binmags+corr3), color='k', s=20)
ax[2].fill_between(time1/60, (binmags+corr3-0.05), (binmags+corr3+0.05), 
                 zorder=0, color='xkcd:light grey')
ax[2].axhline((binmags+corr3).min(), zorder=0, color='xkcd:grey', ls='--')

legend_elements = [Line2D([0], [0], marker='o', color='k', ls='', label='Corrected V\nmagnitude'),
                   Line2D([0], [0], color='xkcd:light grey', lw=10, label='$\sigma(V)$'),
                  Line2D([0], [0], color='xkcd:grey', lw=5, ls='--', label=r'V$_{\rm{max}}$'+'$\simeq${:.3}'.format((binmags+corr3).min()))]
ax[2].legend(handles=legend_elements, loc='lower center', ncol=1, frameon=True, 
          handlelength=2, fontsize=16, markerscale=2)

labels=['TYC 1912-1125-1','GSC 01912-01105','GSC 01912-01015']
for i in range(3):
#     ax[i].set_xlabel('Hours since '+ tobs1[0][11:] + ' UTC', fontsize=20)
    ax[i].set_ylabel('Apparent V magnitude', fontsize=20)
    ax[i].tick_params(labelsize=20)
    ax[i].invert_yaxis()
    ax[i].set_title("Light Curve Corrected with {}".format(labels[i]), fontsize=20)
fig.supxlabel('Hours since '+ tobs1[0][11:] + ' UTC', fontsize=30, y=0.08)    
# plt.subplots_adjust(wspace=0.3)
plt.savefig("Cal_Stars_Light_Curve.jpg")


# In[ ]:





# In[73]:


# fig, ax = plt.subplots(figsize=(10,8))
# mag_corr = binmags+dVavgs
# # maxmag = (binmags+dVavgs).max()
# maxmag = mag_corr.max()

# binmags = - 2.5 * np.log10(binary.aperture_sum_bkgsub_local.values.astype(float))
# plt.scatter(time1/60/24, (binmags+dVavgs), color='k', s=20)
# plt.fill_between(time1/60/24, (binmags+dVavgs-dVhavg), (binmags+dVavgs+dVhavg), 
#                  zorder=0, color='xkcd:light grey')
# plt.xlabel(tobs1[0][:10] + ': ' + 'Days since '+ tobs1[0][11:] + ' UTC', fontsize=20)
# plt.ylabel('Apparent V magnitude', fontsize=20)
# ax.tick_params(labelsize=18)
# plt.axvline(time1[15]/60/24, color='r', lw=1, ls='--')
# plt.fill_between(np.linspace(time1[13],time1[17],10)/60/24, 11.5,12.5, zorder=0, color='r', alpha=0.075)
# plt.axvline(149.7/60/24, color='r', lw=1, ls='--')
# plt.axhline(mag_corr.max(), zorder=0, color='xkcd:grey', ls='--')
# plt.fill_between(np.linspace(time1[60],time1[68],10)/60/24, 11.5,12.5, zorder=0, color='r', alpha=0.075)
# #plt.axvline(0.16, color='r', lw=1, ls='--')
# #plt.fill_between(np.linspace(0.155,0.18,10), 11.5,12.5, zorder=0, color='r', alpha=0.1)
# legend_elements = [Line2D([0], [0], marker='o', color='k', ls='', label='Corrected V\nmagnitude'),
#                    Line2D([0], [0], color='xkcd:light grey', lw=10, label='$\sigma(V)$'),
#                     Line2D([0], [0], color='r', lw=5, ls='--', label='orbit extremas'),
#                   Line2D([0], [0], color='r', lw=10, alpha=0.1, label=r'$\sigma(\rm{extremas})$'),
#                   Line2D([0], [0], color='xkcd:grey', lw=5, ls='--', label=r'V$_{\rm{max}}$'+'$\simeq${:.3}'.format(mag_corr.max()))]
# # plt.ylim(11.5,12.5)
# ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor= (1.01, 0.45), ncol=1,
#             borderaxespad=0, frameon=True, handlelength=3, fontsize=16, markerscale=2)


# In[74]:


time1[63]/60/24


# In[75]:


(time1[67]+time1[62])/2/60/24


# In[76]:


mag_corr = binmags + corrs


# In[77]:


from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

x = time1[0:33]/60 #np.arange(-30, 31, step) #
y = mag_corr[0:33]#azP

n = len(x)                          #the number of data
mean = np.mean(x)            #note this correction
# sigma = sum(y*(x-mean)**2)/n        #note this correction
sigma = np.std(x)

def gaus(x,a,x0,sigma,c):
#     return a * (x-x0)**2 +c 
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + c

popt2,pcov2 = curve_fit(gaus,x,y, p0=[1, mean, sigma, 14.2]) 
print(popt2)

fig, ax = plt.subplots(figsize=(10,8))
plt.scatter(x,y,color='k',label='data',s=10)
plt.plot(np.linspace(x.min(),x.max(),1000),
         gaus(np.linspace(x.min(),x.max(),1000),*popt2),color='r',ls='--',label='fit', lw=1.0)
ax.legend(loc="best", fontsize=20, markerscale=2)
ax.invert_yaxis()
plt.title('Fit for "peak" around t = {:.4} hours'.format(time1[15]/60), fontsize=25)
plt.xlabel('Time (days)', fontsize=20)
plt.ylabel('Apparent V magnitude', fontsize=20)
plt.savefig("day1_1st_lc_peak_fit.jpg")
plt.show()


# In[78]:


pcov2
perr2 = np.sqrt(np.diag(pcov2))
perr2


# In[79]:


from scipy.optimize import curve_fit
from scipy import asarray as ar,exp

x = time1[15:110]/60/24 #np.arange(-30, 31, step) #
y = mag_corr[15:110]#azP

n = len(x)                          #the number of data
mean = np.mean(x)            #note this correction
print(mean)
# sigma = sum(y*(x-mean)**2)/n        #note this correction
sigma = np.std(x)

def gaus(x,a,x0,sigma,c):
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + c

popt,pcov = curve_fit(gaus, x, y, p0=[-1, mean, sigma, 14.6])
print(popt)

fig, ax = plt.subplots(figsize=(10,8))
plt.scatter(x,y,color='k',label='data',s=10)
plt.plot(np.linspace(x[0], x[-1],1000),
         gaus(np.linspace(x[0], x[-1],1000),*popt),color='r',ls='--',label='fit', lw=1.0)
plt.legend(loc="best", fontsize=20, markerscale=2)
ax.invert_yaxis()
plt.title('Fit for peak Centred Near 2.5 Hours', fontsize=25)
plt.xlabel('Time (days)', fontsize=20)
plt.ylabel('Apparent V magnitude', fontsize=20)
plt.savefig("day1_2nd_lc_peak_fit.jpg")
plt.show()


# In[80]:


pcov
perr = np.sqrt(np.abs(np.diag(pcov)))
perr


# In[81]:


(time1[62]-time1[15])/60/24*4


# In[96]:


# Calculating Luminosities of stars from Light curve data
au = 1.495978707e11 # m
p = 2.74e-3
up = 0.02e-3

A_plaskett = np.pi * (1.8 / 2) ** 2
D_bin = 206264 / p * au
uD_bin = (206264 / p) *  up * au

K = A_plaskett / D_bin ** 2
uK = K * 2 / (D_bin) * uD_bin

Vmin1 = mag_corr.max()
uVmin = dC
Vmin2 = mag_corr[0:33].max()
Vmax = mag_corr.min()

L1 = 4 / K * 10 ** (-1 * Vmin1 / 2.5) # W
uL1 = np.sqrt((L1 / K * uK) ** 2 + (4 / K * 2 * np.log(10) / (5 * 10 ** (-1 * Vmin1 / 2.5)) * uVmin) ** 2) # W

L2 = 4 / K * (10 ** (-1 * Vmax / 2.5) - 10 ** (-1 * Vmin1 / 2.5)) # W

uL2 = np.sqrt((L2 / K * uK) ** 2 + (4 / K * 2 * np.log(10) / (5 * 10 ** (-1 * Vmin1 / 2.5)) * uVmin) ** 2 + (4 / K * 2 * np.log(10) / (5 * 10 ** (-1 * Vmax / 2.5)) * uVmin) ** 2) # W

f = open("L_and_M.txt", 'w')
L = []
L1string = "The Luminosity of 1 is {} \pm {} Watts \n".format(L1, uL1)
Vm1string = "The first V min is {} \pm {} \n".format(Vmin1, uVmin)
L2string = "The Luminosity of 2 is {} \pm {} Watts \n".format(L2, uL2)
Vm2string = "The second V min is {} \pm {} \n".format(Vmin2, uVmin)
Vmaxstring = "Max V is {} \pm {} \n".format(Vmax, uVmin)
D_binstring = "The distance to the binary is {} \pm {} metres \n".format(D_bin, uD_bin)
Kstring = "The constant K is {} \pm {} \n".format(K, uK)

L.append(D_binstring)
L.append(Kstring)
L.append(Vm1string)
L.append(Vm2string)
L.append(Vmaxstring)
L.append(L1string)
L.append(L2string)


Lsun = 3.828e26 # W
Msun = 1.9885e30 # kg

M1 = Msun * (L1 / (1.4 * Lsun)) ** (2 / 7)
uM1 =  (2 / 7) * Msun * (L1 / (1.4 * Lsun)) ** (-5 / 7) * (1 / (1.4 * Lsun)) * uL1

M2 = Msun * (L2 / (1.4 * Lsun)) ** (2 / 7)
uM2 =  (2 / 7) * Msun * (L2 / (1.4 * Lsun)) ** (-5 / 7) * (1 / (1.4 * Lsun)) * uL2

q = M2 / M1
uq = np.sqrt((1 / M1 * uM2) ** 2 + (-1 * M2 / (M1 ** 2) * uM2) ** 2)

M1string = "Star 1's mass is approximately {} \pm {} kg \n".format(M1, uM1)
M2string = "Star 2's mass is approximately {} \pm {} kg \n".format(M2, uM2)

qstring = "the ratio of the masses, q = M2 / M1 is {} \pm {} \n".format(q, uq)

L.append(M1string)
L.append(M2string)
L.append(qstring)

f.writelines(L)
f.close()


# In[97]:


print(L)


# ## Day 2:

# In[37]:


aimgf2.summary


# In[38]:


# Initial guess:
target2 = (590,1125)
comp2 = [(332,1166), (589,1406), (636,1523)]
ourstars2 = [(590,1125), (332,1166), (589,1406), (636,1523)]


# In[39]:


exproc = CCDData.read(aimg_path2/aimgf2.files[0], unit='adu')
fig, ax = plt.subplots(figsize=(10,10))
show_image(exproc, cmap='Greys', percl=99.9, fig=fig, ax=ax, show_colorbar=False)
apertures = CircularAperture(target, r=11.);
apertures.plot(color='r', lw=2.5, alpha=1);
apertures = CircularAperture(target, r=12.);
apertures.plot(color='xkcd:black', ls='--', lw=1.5, alpha=1);
apertures = CircularAperture(target, r=19.);
apertures.plot(color='xkcd:black', lw=1.5, ls='--', alpha=1);

rvals = [11, 11, 11, 11]
for r, compstar in enumerate(comp2):
    apertures = CircularAperture(compstar, r=rvals[r]);
    apertures.plot(color='xkcd:bright blue', lw=3, alpha=1);
    apertures = CircularAperture(compstar, r=rvals[r]+2);
    apertures.plot(color='xkcd:bright pink', lw=2, ls='--', alpha=1);
    apertures = CircularAperture(compstar, r=rvals[r]+8);
    apertures.plot(color='xkcd:bright pink', lw=2, ls='--', alpha=1);
#plt.xlim(250,750)
plt.ylim(top=2200)
ax.tick_params(labelsize=18)
plt.xlabel('X pixel', fontsize=20)
plt.ylabel('Y pixel', fontsize=20)

legend_elements = [Line2D([0], [0], lw=3, color='xkcd:bright blue', ls='-', label='Calib. star\naperture'),
                   Line2D([0], [0], lw=3, color='xkcd:bright pink', ls='--', label='Calib. star\nBG annulus'),
                   Line2D([0], [0], color='r', lw=3, label='Target aperture'),
                  Line2D([0], [0], color='k', lw=3, label='Target BG annulus')]

ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor= (1.02, 0.4), ncol=1,
            borderaxespad=0, frameon=True, handlelength=3, fontsize=16)


# In[40]:


# Zoomed in:

fig, ax = plt.subplots(figsize=(10,10))
show_image(exproc, cmap='Greys', percl=99.9, fig=fig, ax=ax, show_colorbar=False)
apertures = CircularAperture(target, r=11.);
apertures.plot(color='r', lw=2.5, alpha=1);
apertures = CircularAperture(target, r=12.);
apertures.plot(color='k', ls='--', lw=2, alpha=1);
apertures = CircularAperture(target, r=19.);
apertures.plot(color='k', lw=2, ls='--', alpha=1);

rvals = [11, 11, 11, 11]
for r, compstar in enumerate(comp):
    apertures = CircularAperture(compstar, r=rvals[r+1]);
    apertures.plot(color='xkcd:bright blue', lw=3, alpha=1);
    apertures = CircularAperture(compstar, r=rvals[r+1]+2);
    apertures.plot(color='xkcd:bright pink', lw=2, ls='--', alpha=1);
    apertures = CircularAperture(compstar, r=rvals[r+1]+8);
    apertures.plot(color='xkcd:bright pink', lw=2, ls='--', alpha=1);
plt.xlim(250,750)
plt.ylim(1000,1650)
ax.tick_params(labelsize=18)
plt.xlabel('X pixel', fontsize=20)
plt.ylabel('Y pixel', fontsize=20)

legend_elements = [Line2D([0], [0], lw=3, color='xkcd:bright blue', ls='-', label='Calib. star\naperture'),
                   Line2D([0], [0], lw=3, color='xkcd:bright pink', ls='--', label='Calib. star\nBG annulus'),
                   Line2D([0], [0], color='r', lw=3, label='Target aperture'),
                  Line2D([0], [0], color='k', lw=3, label='Target BG annulus')]

ax.legend(handles=legend_elements, loc='lower left', bbox_to_anchor= (0.01, 1.01), ncol=2,
            borderaxespad=0, frameon=True, handlelength=3, fontsize=16)


# In[41]:


phot_table


# In[42]:


aimg_path2


# In[43]:


import pandas as pd
data2 = pd.DataFrame(columns=['time', 'star', 'aperture_sum', 'total_bkg_local', 
                             'aperture_sum_bkgsub_local', 'aperture_sum_bkgsub_global', 'aperstats_sum'])
apertures = CircularAperture(ourstars2, r=11.);
annulus_apertures = CircularAnnulus(ourstars2, r_in=12., r_out=19.)

for i in range(len(aimgf2.files)):
    exproc = CCDData.read(aimg_path2/aimgf2.files[i], unit='adu')

    mask = make_source_mask(exproc, nsigma=3, npixels=5, dilate_size=11)
    mean, median, std = sigma_clipped_stats(exproc, sigma=3.0, mask=mask)
    phot_table_local = aperture_photometry(exproc - median, apertures)
    
    aperstats = ApertureStats(exproc.data, annulus_apertures)
    sums = aperstats.sum
    bkg_mean = aperstats.median
    phot_table = aperture_photometry(exproc, apertures)
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output
    aperture_area = apertures.area_overlap(exproc)
    total_bkg = bkg_mean * aperture_area
    phot_bkgsub = phot_table['aperture_sum'].value - total_bkg

    phot_table['total_bkg'] = total_bkg
    phot_table['aperture_sum_bkgsub'] = phot_bkgsub
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output
    
    dat = np.array([time2[i]*np.ones(4), phot_table['id'].value.astype(int), phot_table['aperture_sum'].value.astype(float),
                    phot_table['total_bkg'].value.astype(float), phot_table['aperture_sum_bkgsub'].value.astype(float), 
                    phot_table_local['aperture_sum'].value.astype(float), sums]).T
    data2 = data2.append(pd.DataFrame(dat, columns=data2.columns))


# In[44]:


binary2 = data2.loc[data2.star==1]
comp2 = data2.loc[data2.star.isin([2,3,4])]


# In[45]:


plt.scatter(binary2.time.values, -2.5*np.log10(binary2.aperture_sum_bkgsub_local.astype(float)), s=10)


# In[46]:


mjdt1 = np.array(aimgf1.summary['mjd-obs'].value)
dVavgs = []
fig, ax = plt.subplots(figsize=(10,8))
for i,t in enumerate(time2):
    dtime = comp2.loc[comp2.time.values==t]
    dV = dtime.Vmag.values - (2.5*np.log10(dtime.aperstats_sum.values.astype(float)))
    dVl = dtime.Vmag.values - dtime.uVmag.values - (2.5*np.log10(dtime.aperstats_sum.values.astype(float)))
    dVh = dtime.Vmag.values + dtime.uVmag.values - (2.5*np.log10(dtime.aperstats_sum.values.astype(float)))
    dVavg = np.mean(dV)
    dVavgs.append(dVavg)
    dVhavg = np.mean(dVh) - dVavg
    dVlavg = dVavg - np.mean(dVl)
    binmag = 2.5 * np.log10(binary2.loc[binary2.time.values==t].aperstats_sum.values.astype(float))
    #plt.scatter(t/60, binmag+dVavg, color='k', s=10)
binmags = 2.5 * np.log10(binary2.aperstats_sum.values.astype(float))
dVavgs = np.array(dVavgs)
maxmag = (binmags+dVavgs).max()
plt.scatter(time2/60, (binmags-dVavgs), color='k', s=20)
#plt.fill_between(time2/60, (binmags+dVavgs-dVhavg), (binmags+dVavgs+dVhavg), zorder=0, color='xkcd:light grey')
plt.xlabel(tobs1[0][:10] + ': ' + 'Hours since '+ tobs1[0][11:] + ' UTC', fontsize=20)
plt.ylabel('Apparent V magnitude', fontsize=20)
ax.tick_params(labelsize=18)

legend_elements = [Line2D([0], [0], marker='o', color='k', ls='', label='Corrected V\nmagnitude'),
                   Line2D([0], [0], color='xkcd:light grey', lw=10, label='$\sigma(V)$')]

#ax.legend(handles=legend_elements, loc='upper left', fontsize=18,markerscale=2)


# In[ ]:


comp2['Amag'] = -2.5 * np.log10(comp2.aperstats_sum.values.astype(float));
comp2['Vmag'] = np.concatenate(np.repeat([Vstars], 39, axis=0));
comp2['uVmag'] = np.concatenate(np.repeat([uVstars], 39, axis=0));

