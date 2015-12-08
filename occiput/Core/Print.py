
# occiput  
# Stefano Pedemonte
# Harvard University, Martinos Center for Biomedical Imaging 
# April. 2015, Boston, MA

from math import pi 


def millisec_to_min_sec(ms):
    """Convert millisecond [integer] to string: min sec msec."""
    totsec = int(ms)/1000
    sec = totsec%60
    min = int(ms)/1000/60 
    msec = ms%(1000)
    return str(min).zfill(2)+" [min]  "+ str(sec).zfill(2)+" [sec]  " + str(msec).zfill(3)+" [msec]"

def pretty_print_large_number(number):
    """Given a large number, it returns a string of the sort:  '10.5 Thousand' or '12.3 Billion'. """
    s = str(number).ljust(12)
    if number > 0 and number < 1e3: 
        pass 
    elif number >= 1e3 and number < 1e6: 
        s = s + "  (%3.1f Thousand)"%(number*1.0/1e3)
    elif number >= 1e6 and number < 1e9: 
        s = s + "  (%3.1f Million)"%(number*1.0/1e6)
    elif number >= 1e9 and number < 1e12: 
        s = s + "  (%3.1f Billion)"%(number*1.0/1e9)
    elif number >= 1e12 and number < 1e15: 
        s = s + "  (%3.1f Trillion)"%(number*1.0/1e12)
    return s

def print_percentage(number): 
    """Given a number between 0 and 1, it returns a string indicating percentage: e.g. 10.2 %%"""
    return "%2.2f %%"%((1.0*number)*100) 

def rad_to_deg(rad): 
    """Convert radians to degrees."""
    return rad * 180.0 / pi

def deg_to_rad(deg): 
    """Convert degrees to radians."""
    return deg * pi / 180.0 

def array_to_string(array, format="%3.3f "): 
    """Return a string from an array, with given formatting. """
    s = ""
    for i in array: 
        s = s+format%i 
    return s

    