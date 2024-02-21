# 

def Find_stars_photometery_v4(file, number_of_stars,FWHM_x,saturation_file_name,csv_filename):
    
    from astropy.stats import sigma_clipped_stats
    from photutils.detection import DAOStarFinder
    import numpy as np
    from photutils.aperture import ApertureStats
    from photutils.aperture import CircularAperture, CircularAnnulus
    from photutils.aperture import aperture_photometry
    from astropy.io import fits

        
    threshold_for_DAO = 5.0
    
    hdul = fits.open(file)    
    data = hdul[0].data
    header_data = hdul[0].header
    
    arcsec_pix = header_data['SCALE']
    exposure_time = header_data['EXPTIME']
    
    # If we get a new camera, we need to add it here!
    if header_data['INSTRUME'] == 'ZWO CCD ASI1600MM Pro':
        N_bits = 12
        read_noise = 1.5
    elif header_data['INSTRUME'] == 'ZWO CCD ASI2600MM Pro':
        N_bits = 16
        read_noise = 3.5

    
    # Our usual sigma for a focused star is about 6 arc-secons, need to convert to pixels
    # here the arcsex per pix is relative to the original image, if we've binned it we need to take that into account
    # for the binned image
    
    FWHM = 6 /  arcsec_pix

    mean, median, std = sigma_clipped_stats(data, sigma=3.0)

        # (1) Use the first guesses of FWHM etc. to find "1st draft" stars

    daofind = DAOStarFinder(fwhm=FWHM, threshold=threshold_for_DAO*std, brightest=5,exclude_border=True)  
    sources = daofind(data - median)  
    for col in sources.colnames:  
        sources[col].info.format = '%.8g'  # for consistent table output
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

    print('Done 1st draft find, found '+ str(len(positions))+' stars')

        # (2) Run aperature photometry on the "1st draft" stars
        # Note from VS thesis the aperature is about 5*FWHM to ensure all star is measured

    aperture = CircularAperture(positions, r=FWHM*1.3)
    aperstats = ApertureStats(data - median, aperture) #So that 2D Gaussian is sensible
    new_fwhm = aperstats.fwhm.value
        #print(new_fwhm)
    new_FWHM = np.nanmean(new_fwhm)

    print('Done 1st draft photometry')
    
    
    print('Old FWHM = %.2f pixs' % FWHM)
    FWHM_Arc = FWHM*arcsec_pix
    print('Old FWHM = %.2f arcsec' % FWHM_Arc) 
    
    print('New FWHM = %.2f pixs' % new_FWHM)
    
    FWHM_Arc = new_FWHM*arcsec_pix
    print('New FWHM = %.2f arcsec' % FWHM_Arc)
        # (3) Update the FWHM and re-run star finder

    FWHM = new_FWHM

    daofind = DAOStarFinder(fwhm=FWHM, threshold=threshold_for_DAO*std, brightest=number_of_stars,exclude_border=True)  
    sources = daofind(data - median)  
    for col in sources.colnames:  
        sources[col].info.format = '%.8g'  # for consistent table output
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

    print('Done 2nd round find, found '+ str(len(positions))+' stars')

        # (4) Re-run aperature photometery with the new star list.
        # Set aperature to be 3*FWHM and anulus to have some area.

    r_value     = FWHM*FWHM_x
    # To get an outer annulus with 4x the area of the circular aperture, use this. Some way to ensuring tte uncertainty in the
    # mean or median sky is as low as we can get, as the spread of the sky in the actual aperture is always going to be the limiting facotr.
    # Excepet, if we add in to much annulus this will make our errors larger as we add in more and more read noise. Oh bugger.
    r_out_value = r_value*np.sqrt(5)
    r_in_value = r_value*2

    ANNULUS_aperture = CircularAnnulus(positions, r_in=r_in_value, r_out=r_out_value)
    CIRCULAR_aperture = CircularAperture(positions, r=r_value)

    ANNU_stats  = ApertureStats(data, ANNULUS_aperture)
    CIRC_stats = ApertureStats(data, CIRCULAR_aperture)

    phot_table = aperture_photometry(data, CIRCULAR_aperture)
    phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
    df_photometry = phot_table.to_pandas()


    df_photometry['Circ_sum_flux']      = CIRC_stats.sum
    df_photometry['Circ_area (reduced pixs)']          = CIRC_stats.sum_aper_area.value
    df_photometry['Max_Flux']           = CIRC_stats.max

    df_photometry['FWHM_from_aperature (reduced pix)'] = FWHM
    df_photometry['FWHM_from_aperature (arcsec)'] = FWHM * arcsec_pix

    df_photometry['Annu_sum_flux']           = ANNU_stats.sum
    df_photometry['Annu_area (reduced pix)'] = ANNU_stats.sum_aper_area.value
    df_photometry['Annu_pix_flux_median']    = ANNU_stats.median

    df_photometry['Circ_sum_flux_minus_back'] = \
        df_photometry['Circ_sum_flux'] - df_photometry['Annu_pix_flux_median'] * df_photometry['Circ_area (reduced pixs)']

    df_photometry.sort_values(by=['Circ_sum_flux_minus_back'],ascending=False,inplace=True,ignore_index=True)

    df_photometry['Apperature_radius (reduced pix)'] = r_value
    df_photometry['Anulus_in_radius (reduced pix)'] = r_in_value
    df_photometry['Anulus_out_radius (reduced pix)'] = r_out_value

    df_photometry['Apperature_radius (arcsec)'] = r_value * arcsec_pix
    df_photometry['Anulus_in_radius (arcsec)'] = r_in_value* arcsec_pix
    df_photometry['Anulus_out_radius (arcsec)'] = r_out_value * arcsec_pix
    
    #### SATURATION??
    
   # saturation_file_name = file.replace('reduced.fits', 'saturated.fits')
    saturation_image = fits.getdata(saturation_file_name)
    
    ANNU_stats_saturated  = ApertureStats(saturation_image, ANNULUS_aperture)
    CIRC_stats_saturated = ApertureStats(saturation_image, CIRCULAR_aperture)
   
    is_saturated = ANNU_stats_saturated.min * CIRC_stats_saturated.min
    
    

    df_photometry['Saturated'] = is_saturated
    number_saturated  = (df_photometry['Saturated'] == 0).sum()
    ###########
    print("Number of saturated stars = ",str(number_saturated))
    
    df_photometry['Number_FWHM'] = FWHM_x

    df_photometry['Exposure_sec'] = exposure_time

    # Compute uncertainty on the source flux

    Er1 = df_photometry['Circ_sum_flux']/df_photometry['Exposure_sec']
    
    Er2 = (df_photometry['Circ_area (reduced pixs)'] * df_photometry['Annu_sum_flux']  )/  (df_photometry['Exposure_sec']*df_photometry['Annu_area (reduced pix)']  )

    #Read noise is camera specific!
    Er3 = 2*df_photometry['Circ_area (reduced pixs)']*(header_data['XBINNING']**2)* read_noise**2 / df_photometry['Exposure_sec']**2


    Er = np.sqrt(Er1+Er2+Er3)

    df_photometry['Circ_sum_flux_minus_back_ERROR'] = Er

    mag_error = -2.5 * Er / (df_photometry['Circ_sum_flux_minus_back'] * np.log(10))

    df_photometry['Source_Mag'] = -2.5*np.log10(df_photometry['Circ_sum_flux_minus_back'])
    df_photometry['Source_Mag_Error'] = mag_error

    air_mass = header_data['AIRMASS']
    azimuth =  header_data['OBJCTAZ']
    altitude = header_data['OBJCTALT']
    obs_time     = header_data['DATE-OBS']
    exp_time     = header_data['EXPTIME ']
    filter_colour = header_data['FILTER  ']


    df_photometry['Air mass'] = air_mass
    df_photometry['Azimuth'] = azimuth
    df_photometry['Altitude'] = altitude
    df_photometry['Date Observation'] = obs_time
    df_photometry['Exposure time /s '] = exp_time
    df_photometry['Filter'] = filter_colour
    df_photometry['File'] = file
    df_photometry['BINNED'] = header_data['XBINNING']
    df_photometry['Arcsec_per_pix_binned_image'] = arcsec_pix
    
    #remove stars that are too dim, or negative!
    df_photometry.drop(df_photometry[df_photometry['Circ_sum_flux_minus_back'] <= 100].index, inplace = True)


    print('Done 2nd round photometry')
    print(csv_filename)
    df_photometry.to_csv(csv_filename)
    
    return df_photometry


def create_or_empty_folder(folder_path):
    # Check if the folder exists
    import os
    if os.path.exists(folder_path):
        # If it exists, empty its contents
        for file_name in os.listdir(folder_path):
            file_path = os.path.join(folder_path, file_name)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    os.rmdir(file_path)
            except Exception as e:
                print(f"Error: {e}")
        print(f"Folder '{folder_path}' exists and has been emptied.")
    else:
        # If it doesn't exist, create the folder
        os.makedirs(folder_path)
        print(f"Folder '{folder_path}' has been created.")
        
def update_fits_header(file_path, keyword, new_value, comment="Updated with Astropy"):
    # Open the FITS file
    from astropy.io import fits
    with fits.open(file_path, mode='update') as hdul:
        # Access the header of the primary HDU (HDU = Header Data Unit)
        header = hdul[0].header

        # Update the specified keyword with the new value
        header[keyword] = (new_value, comment)

        # Update the file
        hdul.flush()