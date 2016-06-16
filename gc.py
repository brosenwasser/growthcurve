os.chdir(~/iraf)
import pyfits
from pyraf import iraf
iraf.noao(_doprint=0)
iraf.digiphot(_doprint=0)
iraf.apphot(_doprint=0)
  # what is this filter for?
  #ref_filter = np.loadtxt(target_dir + '/s_extraction/ref_filter.txt', unpack=True,  dtype='str')

  # get image
  # imlist = glob.glob(target_dir + '/img/*_sci.fits')
  # Finds the full path of reference image from the filter
  # ref_image = [image for image in imlist if ref_filter in image][0]
  
#change directory to growthcurve
os.chdir(~/Desktop/growthcurve)
  #for the test run just assign it
ref_image = 'acs_I_030mas_040_sci.fits'
  
  # Set the necessary variables for photometry on the reference image
ref_exptime = pyfits.getheader(ref_image)['EXPTIME']
    #ref_inst = pyfits.getheader(ref_image)['INSTRUME']
    #ref_inst = ref_inst.lower()
    
    #match = (inst_zp == ref_inst) & (filter_zp == ref_filter)
    #ref_zp = zp_zp[match]
    
    
    # Define coordinate and output photometry file paths
    # the file number should be a variable for the full code
gals_coo = 'objects_040.coo'
gals_mag = 'objects_040.mag'

  # Set iraf parameters for reference filter 
iraf.datapars.epadu = ref_exptime
    
iraf.centerpars.calgorithm = 'centroid'
        
iraf.fitskypars.annulus = 40. 
iraf.fitskypars.dannulus = 1. 
    
iraf.photpars.apertures = '1.0,3.0,5.0,7.0,9.0,11.0,13.0,15.0,17.0,19.0,21.0,23.0,25.0,27.0,29.0,31.0,33.0,35.0,37.0,39.0'
    #iraf.photpars.zmag = ref_zp
    
    
    # run phot for the image
iraf.phot(ref_image, gals_coo,gals_mag)
