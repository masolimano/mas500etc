import numpy as np
import hjson
from astropy import units as u
from astropy.constants import h as PLANCK_H, c as LIGHTSPEED

FNU_UNIT = u.erg / u.s / u.cm ** 2 / u.Hz
FLAM_UNIT = u.erg / u.s / u.cm ** 2 / u.AA
class ETC:
    """
    Exposure time calculator for the MAS500 telescope
    """
    # Load instrument data
    with open('telescope_data/mas500_datasheet.hjson', 'r') as hj:
        info_tel = hjson.load(hj)
        info_tel['diameter'] *= u.cm
        collecting_area = np.pi * (info_tel['diameter'] / 2) ** 2
        info_tel['pixelsize'] *= u.um / u.pixel
        info_tel['pixelscale'] *= u.arcsec / u.pixel
        info_tel['darkcurrent'] *= u.electron / u.pixel ** 2 / u.s
        for axis in info_tel['npixels']:
            info_tel['npixels'][axis] *= u.pixel

        for mode in info_tel['ron']:
            info_tel['ron'][mode] *= u.electron / u.pixel

        for flt in info_tel['filters']:
            info_tel['filters'][flt]['quantum_efficency'] *= u.electron # per photon
            info_tel['filters'][flt]['effective_wavelength'] *= u.AA
            info_tel['filters'][flt]['bandpass'] *= u.AA


    def __init__(self, input_file):
        """
        Parse the input file
        """
        with open(input_file, 'r') as hj:
            self.input_dict = hjson.load(hj)

        if self.input_dict['mode'] == 'time':
            # Validate and assign units
            assert 'filter' in self.input_dict
            assert 'read_out_mode' in self.input_dict
            assert 'binning' in self.input_dict
            assert 'source_type' in self.input_dict
            assert 'magnitude' in self.input_dict
            assert 'sky_brightness' in self.input_dict
            assert 'seeing' in self.input_dict

            self.input_dict['magnitude'] *= u.STmag

            self.input_dict['seeing'] *= u.arcsec

            if self.input_dict['binning'] == '1x1':
                self.bin2 = 1
            elif self.input_dict['binning'] == '2x2':
                self.bin2 = 4
            else:
                raise SyntaxError('Please enter "2x2" or "1x1" for the binning')

            if self.input_dict['source_type'] == 'point':

                # FWHM area of the point source in native pixels^2
                self.source_area = np.pi * (0.5 * self.input_dict['seeing'] / self.info_tel['pixelscale']) ** 2

            else:
                raise NotImplementedError(f"Sources of the type {self.input_dict['source_type']} are not implemented.")


        elif self.input_dict['mode'] == 'noise':
            # Validate
            pass

        else:
            # Raise some error
            pass

    def _compute_source_electron_rate(self):
        """
        Obtain the rate of photoelectrons produced by the point source
        """
        flux = self.input_dict['magnitude'].to(FLAM_UNIT)
        specific_photon_rate = flux * self.info_tel['filters'][self.input_dict['filter']]['effective_wavelength'] / (PLANCK_H * LIGHTSPEED)
        bandpass = self.info_tel['filters'][self.input_dict['filter']]['bandpass']
        efficency = self.info_tel['filters'][self.input_dict['filter']]['quantum_efficency']
        self.source_electron_rate = specific_photon_rate * self.collecting_area * bandpass * efficency
        self.source_electron_rate = self.source_electron_rate.decompose()

    def _compute_sky_electron_rate(self):
        """
        Obtain the rate of photoelectrons produced by the sky *per pixel^2*
        """
        chip_npixels = self.info_tel['npixels']['x'] * self.info_tel['npixels']['y']
        solid_angle = chip_npixels * self.info_tel['pixelscale'] ** 2
        solid_angle = solid_angle.to(u.arcsec ** 2)
        magnitude_in_fov = (self.input_dict['sky_brightness'] - 2.5 * np.log10(solid_angle.value)) * u.STmag
        sky_photon_rate = (magnitude_in_fov.to(FLAM_UNIT) * self.info_tel['filters'][self.input_dict['filter']]['effective_wavelength'] / (PLANCK_H * LIGHTSPEED)).decompose()
        bandpass = self.info_tel['filters'][self.input_dict['filter']]['bandpass']
        efficency = self.info_tel['filters'][self.input_dict['filter']]['quantum_efficency']
        self.sky_electron_rate = sky_photon_rate * self.collecting_area * bandpass * efficency  / chip_npixels


    def _solve_noise_mode(self):
        """
        Enter a exposure time and obtain noise level
        """
        pass

    def _solve_time_mode(self):
        """
        Enter a magnitude and the target SNR and solve for time
        """
        self._compute_source_electron_rate()
        self._compute_sky_electron_rate()
        a = (self.source_electron_rate ** 2).value
        b = -1.0 * self.input_dict['target_snr'] ** 2 * (self.source_electron_rate + \
                                                         self.sky_electron_rate * self.source_area + \
                                                         self.info_tel['darkcurrent'] * self.source_area)
        b = b.decompose().value
        c = -1.0 * self.input_dict['target_snr'] ** 2 * (self.info_tel['ron'][self.input_dict['read_out_mode']]  ** 2)\
                 * self.source_area / self.bin2
        c = c.decompose().value

        exposure_time = (-b + np.sqrt(b ** 2 - 4 * a * c) / (2 * a)) * u.s
        return exposure_time


if __name__ == '__main__':
    etc = ETC('test_input/input_tmode_00.hjson')
    t = etc._solve_time_mode()
    print(t)

