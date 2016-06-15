    # Preparation for running growth curve photometry on reference filter
    ref_filter = np.loadtxt(target_dir + '/s_extraction/ref_filter.txt', unpack=True,  dtype='str')
    # Remove string from array
    ref_filter = ref_filter[()]
    
    inst_zp, filter_zp, zp_zp = np.loadtxt(pydir + '/data/legus_zeropoints.tab', unpack=True, dtype='str')
    
    imlist = glob.glob(target_dir + '/img/*_sci.fits')
    
    # Finds the full path of reference image from the filter
    ref_image = [image for image in imlist if ref_filter in image][0]
    
    # Set the necessary variables for photometry on the reference image
    ref_exptime = pyfits.getheader(ref_image)['EXPTIME']
    ref_inst = pyfits.getheader(ref_image)['INSTRUME']
    ref_inst = ref_inst.lower()
    
    match = (inst_zp == ref_inst) & (filter_zp == ref_filter)
    ref_zp = zp_zp[match]
    
    # ref_zp is a string within an array, so need to turn into a float
    ref_zp = float(ref_zp[0])
    
    # Define coordinate and output photometry file paths
    stars_coo = target_dir + '/photometry/isolated_stars.coo'
    stars_mag = target_dir + '/photometry/isolated_stars.mag'
    
    clusters_coo = target_dir + '/photometry/isolated_clusters.coo'
    clusters_mag = target_dir + '/photometry/isolated_clusters.mag'
     
    # Remove phot files if they already exist
    if os.path.exists(stars_mag) == True:
        os.remove(stars_mag)
           
    if os.path.exists(clusters_mag) == True:
        os.remove(clusters_mag)
    
           
    # Set iraf parameters for reference filter 
    iraf.datapars.epadu = ref_exptime
    
    iraf.centerpars.calgorithm = 'centroid'
        
    iraf.fitskypars.annulus = 21. 
    iraf.fitskypars.dannulu = 1. 
    
    iraf.photpars.apertures = '1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0'
    iraf.photpars.zmag = ref_zp
    
    # For isolated stars
    iraf.phot(ref_image, stars_coo, stars_mag)
    
    # For isolated clusters
    iraf.phot(ref_image, clusters_coo, clusters_mag)
         
    # Move phot results to new files, remove INDEFs
    stars_growth = target_dir + '/photometry/isolated_stars_growth.mag'
    clusters_growth = target_dir + '/photometry/isolated_clusters_growth.mag'
    
    cmd = 'grep "*" ' + stars_mag + ' > ' + stars_growth
    os.system(cmd)
        
    cmd = 'sed -i.bak "s/INDEF/99.999/g" ' + stars_growth
    os.system(cmd)
        
    cmd = 'grep "*" ' + clusters_mag + ' > ' + clusters_growth
    os.system(cmd)

    cmd = 'sed -i.bak "s/INDEF/99.999/g" ' + clusters_growth 
    os.system(cmd)

    # Remove .bak files to prevent confusion
    bak_stars = stars_growth + '.bak'
    bak_clusters = clusters_growth + '.bak'
    os.remove(bak_stars)
    os.remove(bak_clusters)
    
    #####  Growth curves of isolated stars and clusters  #####
          
    # Stars
    aper_st, flux_st = np.loadtxt(stars_growth, unpack=True, usecols=(0,3))
    
    ratio_st = np.empty(len(aper_st))
    
    naper = 20
    nstar = len(aper_st)/naper
    aper_ind = naper - 1
    
    for k in range(nstar):
    
        for i in range(naper):
    
            ratio_st[i + k*naper] = flux_st[i + k*naper]/flux_st[aper_ind + k*naper]
            
    
    # Plot growth curves   
    fig = plt.figure(figsize = (7,7))
       
    aper_x = np.arange(naper) + 1
        
    for i in range(nstar):
        
        ratio_y = ratio_st[i*naper:(i + 1)*naper]
        plt.plot(aper_x, ratio_y, 'y-')
        plt.annotate(str(i + 1), xy=(8.0, ratio_y[7]),
            horizontalalignment='left', verticalalignment='top', fontsize=6)

        
    plt.plot(aper_x, med_st, 'r-' , linewidth=4.0)
    plt.hlines(0.5, 0, 20, color='black', linewidth=2, zorder=10)
        
    plt.ylabel('Normalized Flux ' + ref_filter.upper())
    plt.xlabel('Radius (pix)')
    plt.xlim(1,20)
    plt.minorticks_on()
        
    fig.savefig(target_dir + '/plots/plot_growth_curve.pdf')
