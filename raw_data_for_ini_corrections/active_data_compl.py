import requests
from bs4 import BeautifulSoup as bs
import js2xml
import datetime
import json
# complementary data for 16 countries
# get all information for active cases in https://www.worldometers.info/coronavirus/country/serbia/

def worldmeter_date_to_owid(worldmeter_date):
    datetime_object = datetime.datetime.strptime(worldmeter_date[:3], "%b")
    month_number = str(datetime_object.month) if datetime_object.month>=10 else '0'+str(datetime_object.month)
    owid_date = worldmeter_date[-4:] + '-' + month_number + '-' + worldmeter_date[4:6]
    return owid_date

def complement_active_data(country_name, ISO):
    url = 'https://www.worldometers.info/coronavirus/country/'+country_name+'/'
    soup = bs(requests.get(url).content, 'html.parser')
    active_cases = {}
    for script in soup.find_all('script'):
        if 'Highcharts' in script.text:
            parser = js2xml.parse(script.text)
            date = parser.xpath('//property[@name="categories"]//string/text()')
            title = parser.xpath('//property[@name="title"]//string/text()')
            data = parser.xpath('//property[@name="data"]//number/@value')
            if 'Active Cases' in title:
                for i in range(len(date)):
                        d1, d2 = date[i], data[i]
                        active_cases[worldmeter_date_to_owid(d1)] = d2
    with open(ISO+'_active.json', 'w') as f:
        json.dump(active_cases, f)

# # ---- correct active data for several countries --------------------------------
import time
import random
need_correct_name_list = ['thailand', 'cameroon', 'serbia',
                          'us', 'belgium', 'sweden', 'netherlands', 'ireland', 'uk', 'saint-kitts-and-nevis',
                          'greece', 'spain', 'france', 'finland', 'norway', 'switzerland']
need_correct_ISO_list = ['THA', 'CMR', 'SRB',
                         'USA', 'BEL', 'SWE', 'NLD', 'IRL', 'GBR', 'KNA',
                         'GRC', 'ESP', 'FRA', 'FIN', 'NOR', 'CHE']

for i in range(len(need_correct_name_list)):
    print(i)
    time.sleep(20 + 10 * random.random())
    complement_active_data(need_correct_name_list[i], need_correct_ISO_list[i])
# # -------------------------------------------------------------------------------
