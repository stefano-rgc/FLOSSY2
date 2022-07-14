from astropy import units as u

# ---------------------------------------------------------
# USER INPUT GOES ONLY IN THE SECTIONS BETWEEN LONG HYPHENS
# ---------------------------------------------------------
# Set TIC number of your TESS star
TIC = 349832567
# ---------------------------------------------------------
# Set path containing the prewhitened frequencies. Must be a CSV file with header.
pw_file = f'pw/pw_tic{TIC}.csv'
# Set the only 3 column names used by the code: 'frequency' and 'amplitude'
pw_col_names = {'frequency': 'frequency',
                'frequency_error': 'e_frequency',
                'amplitude': 'amp'}
# Set unis of the 3 columns. Use astropy units.
# See examples at: https://docs.astropy.org/en/stable/units/index.html#module-astropy.units.si
pw_units = {'frequency': u.Hz, # Hz
            'frequency_error': u.Hz, # Hz
            'amplitude': 1e-3*u.dimensionless_unscaled} # ppt
# ---------------------------------------------------------
# Set path containing the periodogram. Must be a CSV file with header.
pg_file = f'pg/pg_tic{TIC}.csv'
# Set the only 2 column names used by the code: 'frequency' and 'amplitude'
pg_col_names = {'frequency': 'freq',
                'amplitude': 'amp'}
# Set unis of the 2 columns. Use astropy units.
# See examples at: https://docs.astropy.org/en/stable/units/index.html#module-astropy.units.si
pg_units = {'frequency': u.Hz, # Hz
            'amplitude': 1e-3*u.dimensionless_unscaled} # ppt
# ---------------------------------------------------------

class PW:
    """Collect prewhitened info"""
    file = pw_file
    col_names = pw_col_names
    units = pw_units
    
class PG:
    """Collect periodogram info"""
    file = pg_file
    col_names = pg_col_names
    units = pg_units
    
class UserInput:
    """Class where to store user inputs"""
    tic = TIC
    pw = PW()
    pg = PG()