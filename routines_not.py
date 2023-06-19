from   astropy import time
from   astropy.coordinates import AltAz, EarthLocation, SkyCoord
import astropy.units as u
from   astropy.coordinates import get_sun, get_moon
from   astroplan import moon
from   misc import bcolors
from   plotsettings import *
from   standard_libraries import *

# Optical
obs_GemminiSouth  = EarthLocation.of_site('gemini_south')
obs_LaSilla       = EarthLocation.of_site('lasilla')
obs_LBT           = EarthLocation.of_site('lbt')
obs_NOT           = EarthLocation.of_site('lapalma')
obs_Palomar       = EarthLocation.of_site('Palomar')
obs_VLT           = EarthLocation.of_site('paranal')

# Radio
obs_ALMA          = EarthLocation.of_site('ALMA')
#obs_ATCA          = astroplan.Observer.at_site('atca')
obs_VLA           = EarthLocation.of_site('vla')

def check_moon(coords, obstime, observatory, avoid=30.*u.degree):

    if not isinstance(avoid, u.Quantity):
        avoid = float(avoid)*u.degree
    else:
        avoid = avoid.to(u.degree)

    _moon = get_moon(obstime, obs_NOT)

    target = coords
    
    moon_alt = _moon.transform_to(AltAz(obstime=obstime, location=observatory)).alt.to(u.deg)
    
    if moon_alt < 0*u.degree:
        print('Moon is down')
        return {'sep': 0, 'fli':np.nan}
    
    else:
        sep = _moon.separation(target)
        fli = moon.moon_illumination(obstime)
        
        msg_sep = 'Moon separation = {sep:.0f} deg, fli = {fli:.0f}%'.format(sep=sep.to(u.degree).value, fli=fli*100)
        #msg_fli = 'Full moon illumination {:.0f}%'.format(fli*100)
        msg_fail = 'Object too close to the moon!'
        
        if sep > avoid:
            print(bcolors.OKGREEN + msg_sep + bcolors.ENDC)
            #print(bcolors.OKGREEN + msg_fli + bcolors.ENDC)
            #print(bcolors.OKGREEN + msg_fail + bcolors.ENDC)
        else:
            print(bcolors.FAIL + bcolors.BOLD + msg_sep + bcolors.ENDC)
            #print(bcolors.FAIL + bcolors.BOLD + msg_fli + bcolors.ENDC)
        
        return {'sep': sep.to(u.degree).value, 'fli':fli*100}

def plot_objects(DATE, CAT, OBS=obs_NOT, SAVE=False):

    # Time array

    midnight_utc = time.Time(DATE, format='isot', scale='utc')+1
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour

    # Plot

    plt.figure(figsize=(9*np.sqrt(2), 9))
    ax = plt.subplot(111)

    frame_time  = AltAz(obstime=midnight_utc + delta_midnight, location=obs_NOT)

    for ii in range(len(CAT)):

        print(bcolors.BOLD + bcolors.HEADER + CAT['NAME'][ii] + bcolors.ENDC)
        
        coords      = SkyCoord(CAT['RA'][ii], CAT['DEC'][ii], unit=(u.hour, u.deg))
        obj_altazs  = coords.transform_to(frame_time)
        obj_airmass = obj_altazs.secz

        moon_distance = check_moon(coords, midnight_utc, obs_NOT)        
        
        mask = (obj_airmass < 5) & (obj_airmass > 0)
        ax.plot(delta_midnight, obj_altazs.alt,
                label='{obj} (moon distance: {sep:.0f}°)'.format(obj=CAT['NAME'][ii], sep=moon_distance['sep']), lw=4)

    # Night

    sunaltazs = get_sun(midnight_utc + delta_midnight).transform_to(frame_time)

    ax.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                    sunaltazs.alt < -0*u.deg, color='navy', zorder=0, alpha=0.2)

    ax.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                    sunaltazs.alt < -18*u.deg, color='navy', zorder=0)

    # Moon

    moonaltazs = get_moon(midnight_utc + delta_midnight).transform_to(frame_time)
    ax.plot(delta_midnight, moonaltazs.alt, lw=2, color='white')

    # Sun

    sunaltazs = get_sun(midnight_utc + delta_midnight).transform_to(frame_time)
    ax.plot(delta_midnight, sunaltazs.alt, lw=2, color='black')


    # Airmass limits

    ax.axhline(90 - np.arccos(1/2.)*180/np.pi, color=vigit_color_12)
    ax.axhline(90 - np.arccos(1/3.)*180/np.pi, color=vigit_color_12, ls='--')

    # Midnight

    ax.axvline(0., color='white', lw=2)

    # Legend

    ax.legend(loc='upper left', fontsize=legend_size-4, bbox_to_anchor=(1, 0.25, 0.5, 0.5))

    # Prettify

    xmin, xmax = -8, 10

    ax.set_xlim(xmin*u.hour, xmax*u.hour)
    labels = ['{:.0f} h'.format(x + 24) if x < 0 else'{:.0f} h'.format(x) for x in (np.arange(xmax-xmin)+xmin)]
    ax.set_xticks((np.arange(18)-8)*u.hour, labels=labels, rotation=45)

    ax.set_ylim(10*u.deg, 90*u.deg)

    ax.set_xlabel('Universal time', fontsize=label_size)
    ax.set_ylabel('Altitude', fontsize=label_size)

    try:
        ax.set_title(DATE + ' → ' + (time.Time(DATE)+1).isot.split('T')[0] + ' (FLI = ' + str(int(moon_distance['fli'])) + '%)', fontsize=label_size)
    except:
        ax.set_title(DATE + ' → ' + (time.Time(DATE)+1).isot.split('T')[0], fontsize=label_size)

    if SAVE:
        plt.savefig('visibility_targets_{}.pdf'.format(DATE))

def plot_standards(DATE, OBS=obs_NOT, SAVE=False):

    # ESO/HST standard stars visible from La Palma
    # http://www.not.iac.es/instruments/alfosc/fluxstandard.html

    list_std       = table.Table(names=('NAME', 'RA', 'DEC'), dtype=('S100', 'S100', 'S100'))
    list_std.add_row(['SP1036+433', '10:39:36.7358', '+43:06:09.212'])
    list_std.add_row(['SP1045+378', '10:48:23.5113', '+37:34:13.083'])
    list_std.add_row(['SP1446+259', '14:49:02.3593', '+25:42:09.192'])
    list_std.add_row(['SP1550+330', '15:51:59.8854', '+32:56:54.329'])
    list_std.add_row(['SP2209+178', '22:11:31.3756', '+18:05:34.177'])
    list_std.add_row(['SP2317-054', '23:19:58.3996', '-05:09:56.171'])

    # Time array

    midnight_utc = time.Time(DATE, format='isot', scale='utc')+1
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour

    # Plot

    plt.figure(figsize=(9*np.sqrt(2), 9))
    ax = plt.subplot(111)

    frame_time  = AltAz(obstime=midnight_utc + delta_midnight, location=OBS)

    for ii in range(len(list_std)):

        print(bcolors.BOLD + bcolors.HEADER + list_std['NAME'][ii] + bcolors.ENDC)
        
        coords      = SkyCoord(list_std['RA'][ii], list_std['DEC'][ii], unit=(u.hour, u.deg))
        obj_altazs  = coords.transform_to(frame_time)
        obj_airmass = obj_altazs.secz
    
        moon_distance = check_moon(coords, midnight_utc, obs_NOT)        

        mask        = (obj_airmass > 0) & (obj_airmass < 5)
        ax.plot(delta_midnight, obj_altazs.alt,
                label='{obj} (moon distance: {sep:.0f}°)'.format(obj=list_std['NAME'][ii], sep=moon_distance['sep']),
                lw=4,
                color=colors_vigit[2*ii])

    # Night

    sunaltazs = get_sun(midnight_utc + delta_midnight).transform_to(frame_time)

    ax.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                    sunaltazs.alt < -0*u.deg, color='navy', zorder=0, alpha=0.2)

    ax.fill_between(delta_midnight, 0*u.deg, 90*u.deg,
                    sunaltazs.alt < -18*u.deg, color='navy', zorder=0)

    # Moon

    moonaltazs = get_moon(midnight_utc + delta_midnight).transform_to(frame_time)
    ax.plot(delta_midnight, moonaltazs.alt, lw=2, color='white')

    # Sun

    sunaltazs = get_sun(midnight_utc + delta_midnight).transform_to(frame_time)
    ax.plot(delta_midnight, sunaltazs.alt, lw=2, color='black')


    # Airmass limits

    ax.axhline(90 - np.arccos(1/2.)*180/np.pi, color=vigit_color_12)
    ax.axhline(90 - np.arccos(1/3.)*180/np.pi, color=vigit_color_12, ls='--')

    # Midnight

    ax.axvline(0., color='white', lw=2)

    # Legend

    ax.legend(loc='upper left', fontsize=legend_size-4, bbox_to_anchor=(1, 0.25, 0.5, 0.5))

    # Prettify

    xmin, xmax = -8, 10

    ax.set_xlim(xmin*u.hour, xmax*u.hour)
    labels = ['{:.0f} h'.format(x + 24) if x < 0 else'{:.0f} h'.format(x) for x in (np.arange(xmax-xmin)+xmin)]
    ax.set_xticks((np.arange(18)-8)*u.hour, labels=labels, rotation=45)

    ax.set_ylim(10*u.deg, 90*u.deg)

    ax.set_xlabel('Universal time', fontsize=label_size)
    ax.set_ylabel('Altitude', fontsize=label_size)

    ax.set_title(DATE + ' → ' + (time.Time(DATE)+1).isot.split('T')[0] + ' (FLI = ' + str(int(moon_distance['fli'])) + '%)', fontsize=label_size)

    if SAVE:
        plt.savefig('visibility_standards_{}.pdf'.format(DATE))
