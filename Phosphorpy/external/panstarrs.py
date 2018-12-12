from astropy.units.quantity import Quantity
from astropy import units as u
from threading import Thread
import urllib
import os

PS_URL = 'http://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos=' \
         '{}+%09{}&filter=color&filter=g&filter=r&filter=i&filter=z&filter=y' \
         '&filetypes=stack&auxiliary=data&auxiliary=mask&size={}' \
         '&output_size=0&verbose=0&autoscale=99.500000&catlist='


def get_all_image_urls(ra, dec, size):
    """
    Extract the image urls from the PAN-STARRS image web-site

    :param ra: The RA component of the coordinate of the target
    :type ra: float
    :param dec: The DEC component of the coordinate of the target
    :type dec: float
    :param size: The size of the image in degree
    :type size: float
    :return: Dict with where the keys are the name of the bands and the items are the urls to images
    :rtype: dict
    """
    bands = ['g', 'r', 'i', 'z', 'y']
    user_agent = ''.join(['Mozilla/5.0 (Macintosh; U; Intel Mac OS X 10_6_4; en-US) ',
                          'AppleWebKit/534.3 (KHTML, like Gecko) Chrome/6.0.472.63 Safari/534.3'])
    headers = {'User-Agent': user_agent}

    if type(size) == Quantity:
        size = size.to(u.arcsec).value/0.25
    else:
        size /= 0.25
    size = int(size)
    # create the url of the panstarrs image website
    url = PS_URL.format(ra, dec, size)

    # call the website and download the source code
    req = urllib.request.Request(url, None, headers)
    response = urllib.request.urlopen(req)
    page = response.read()
    response.close()
    page = str(page)

    # extract the url of the cutout from the source of the image website
    page = page.split('\\n')
    out = {}
    for band in bands:
        img_url = ''
        bband = ''.join(['.', band, '.'])
        for n in page:
            if 'FITS-cutout' in n and bband in n:
                img_url = n
                break
        img_url = img_url.split('Download FITS cutout')[-1]
        img_url = img_url.split('="')[-1]
        img_url = img_url.split('">')[0]
        out.update({band: img_url})
    return out


def _download_band(url, img_path):
    try:
        urllib.request.urlretrieve(url, img_path)
    except urllib.error.URLError:
        pass


def download_all_bands(ra, dec, size, save_path):
    """
    Downloads all five filter-images

    :param ra: The ra component of the coordinate in deg
    :type ra: float
    :param dec: The dec component of the coordinate in deg
    :type dec: float
    :param size: The size of the image in degrees
    :type size: float
    :param save_path: The path to the save place with two parameters
    :type save_path: str
    :returns: Dictionary with the band names as keys and the paths to the downloaded images as items
    :rtype: dict
    """
    # check if the directory exists
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # get the image urls
    img_urls = get_all_image_urls(ra, dec, size)
    ths = []

    path = os.path.join(save_path, 'temp_image_{}.fits')

    out = {}
    # go through all available bands
    for b in img_urls.keys():
        img_path = path.format(b)
        # start a new thread to download every bands
        th = Thread(target=_download_band,
                    args=('http:' + img_urls[b],
                          img_path,))
        th.start()
        ths.append(th)
        out[b] = img_path

    # wait until all bands are finished or until a timeout of 30sec is reached
    for t in ths:
        t.join(timeout=30)

    return out
