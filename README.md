# Exposure Time Calculator for the MAS500 CCD imager at Observatorio El Sauce

Get the code:
```sh
git clone https://www.github/masolimano/mas500etc.git
cd mas500etc
```

Currently the ETC supports two modes. A `time` mode which solves for the exposure time in the CCD equation given a target SNR. And a `noise`  mode which computes the 3-sigma surface brightness sensitivty for a given exposure time.
## Time mode
```py
from etc import ETC

etc = ETC(
    mode='time',
    filter='V',
    source_type='point',
    read_out_mode='slow',
    binning='2x2',
    magnitude=20,  # ST magnitudes
    sky_brightness=21, # ST magnitude per arcsec^2
    target_snr=60,
    seeing=1.0,  # arcsec
)
etc.run()
```

This particular input should output:
```
Exposure time = 575.3 s
```


## Noise mode
```py
etc = ETC(
    mode='noise',
    filter='g',
    source_type='point',
    read_out_mode='slow',
    binning='2x2',
    sky_brightness=21,
    exposure_time=600 # seconds
)
etc.run()
```

This particular input should output:
```
3-sigma surface brightness sensitivty = 8.9e-19 erg / (Angstrom arcsec2 cm2 s)
```

## Read from input file
Additionally, one can provide an input hjson file with the required keywords. Below we show an example input file:
```
{
# --Calculator settings--
  mode: time
  filter: g
  read_out_mode: slow
  binning: 2x2
  target_snr: 15

# --Source properties--
  source_type: point
  magnitude: 20

# --Sky properties--
  sky_brightness: 20
  # in magnitudes / arcsec^2
  seeing: 1.0 
  # in arcseconds
}
```

To run with the calculator with these parameters we use the `read_input` method:
````py
etc = ETC.read_input(input_filename)
etc.run()
```
