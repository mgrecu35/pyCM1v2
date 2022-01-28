url='http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&YEAR=1999&MONTH=07&FROM=0800&TO=0800&STNM=72357'

import urllib.request
import requests

#fp = urllib.request.urlopen(url)
#mybytes = fp.read()
#mystr = mybytes.decode("utf8")
#fp.close()

#print(mystr)


url_get = requests.get(url)
htmltext = url_get.text
#print(htmltext)

from bs4 import BeautifulSoup
soup = BeautifulSoup(htmltext, 'lxml')
