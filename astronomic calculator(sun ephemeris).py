from math import *
import os
import datetime
import pysolar.solar as ps
from geopy.geocoders import Nominatim

#import pytz

# Banner
while True:
    current_time = datetime.datetime.now()
    def my_banner():
            print("+================================================================================================================+")
            print("-----------------------------------------Astronomical Calculator--------------------------------------------------")
            print("---------------------------------------------Sun Ephemeris--------------------------------------------------------")
            print("-----------------------------------Created by : AHMAD FAIZ BIN IBRAHIM--------------------------------------------")
            print("------------------------------------------On 18 February 2023-----------------------------------------------------")
            print("---------------------------------------" + str(current_time) + "--------------------------------------------------")
            print("+================================================================================================================+\n\n")
            return

    my_banner()
    # Input from Observer 
    Latitude = float(input('Latitude: '))
    Longitude = float(input('Longitude: '))

    geolocator = Nominatim(user_agent="my-application/1.0")
    location = geolocator.reverse(f"{Latitude}, {Longitude}")

    Time_Zone = int(input('GMT: '))
    DST_Hour = int(input('DST: '))

    Year = int(input('Year: '))
    Month = int(input('Month: '))
    Day = int(input('Day: '))
    Hour = int(input('Hour: '))
    Minute = int(input('Minute: '))
    Second = int(input('Second: '))
    temperature = int(input('Temperature in Celcius: '))
    ketinggian_dari_aras_laut = int(input('Height above sea level in meters: '))
    humidityinp = float(int(input('Humidity in percentage (For example: 77): ')) / 100)

    # Time related Parameters, Greenwich Mean Sidereal Time & the Sun's Mean Longitude
    ST_hrs = float(Hour + (Minute/60) + (Second/3600))
    UTC_Uncorrected = ST_hrs - Time_Zone - DST_Hour
    UTC_Corrected = UTC_Uncorrected % 24
    UTC_Corrected_d = UTC_Corrected * 15
    UTC_Corrected_r = radians(UTC_Corrected_d)

    def temporary_value():
        a = 0
        if (UTC_Uncorrected < 0):
            a += -1
        
        elif (UTC_Uncorrected > 24):
            a += 1
        
        else:
            a += 0
        
        aaa = a
        bbb = (367 * Year) - 730531.5
        ccc = int((7. * int(Year+ (Month+9)/12))/4)
        ddd = int(275 * Month/9) + Day
        J_200 = aaa+bbb-ccc+ddd
        return J_200

    Days_since_Midnight = UTC_Corrected / 24
    Days_to_0000am_since_Epoch = temporary_value()
    Days_to_Now_since_Epoch = Days_to_0000am_since_Epoch + Days_since_Midnight
    Julian_Centuries_2000 = Days_to_Now_since_Epoch / 36525

    Greenwitch_Mean_Sideral_Time = (6.697374558 +0.06570982441908 *Days_to_0000am_since_Epoch +1.00273790935 * UTC_Corrected +0.000026 * Julian_Centuries_2000)% 24
    Greanwitch_Mean_Sideral_Time_d = Greenwitch_Mean_Sideral_Time * 15

    Sun_Mean_Longitude_d = (Greanwitch_Mean_Sideral_Time_d - 180 -UTC_Corrected_d) % 360
    Sun_Mean_Longitude_r = radians(Sun_Mean_Longitude_d)


    # Astronomical Facts
    Perihelion_Longitude_d = 248.54536+0.017196*Year
    Perihelion_Longitude_r = radians(Perihelion_Longitude_d)

    Ecentricity = 0.017585 - 0.438 * Year/1000000

    Obliquity_d = 23.6993 - 0.00013 * Year
    Obliquity_r =radians(Obliquity_d)

    # Solving Kepler's Theorem & Sun's True Longitude 

    Mean_Anomaly_r = Sun_Mean_Longitude_r - Perihelion_Longitude_r
    Ecentric_Anomaly_r = Mean_Anomaly_r - sin(Mean_Anomaly_r)/(cos(Mean_Anomaly_r)-1/Ecentricity)
    True_Anomaly = 2 * atan(tan(Ecentric_Anomaly_r/2)*sqrt((1+Ecentricity)/(1-Ecentricity)))
    Sun_True_Longitude_r = True_Anomaly + Perihelion_Longitude_r
    Sun_True_Longitude_d = degrees(Sun_True_Longitude_r)

    # Sun's Declination, Right Ascension & the Equation of Time 

    Sun_Declination_r = asin(sin(Obliquity_r)*sin(Sun_True_Longitude_r))
    Sun_Declination_d = degrees(Sun_Declination_r)

    Right_Ascension_r = atan2(cos(Obliquity_r)*sin(Sun_True_Longitude_r),cos(Sun_True_Longitude_r))
    Right_Ascension_d = degrees(Right_Ascension_r) % 360
    Right_Ascension = Right_Ascension_d / 15

    Equation_of_Time_d = Greanwitch_Mean_Sideral_Time_d - Right_Ascension_d - UTC_Corrected_d + 180

    def EoTas():
        EoTx = Equation_of_Time_d
        if (EoTx<-180):
            EoTx+=360

        elif (EoTx>180):
            EoTx -= 360
            
        else:
            EoTx+=0

        return EoTx

    Equation_of_Time_Astro_d = EoTas()
    Equation_of_Time_Gnomical_d = - Equation_of_Time_Astro_d
    Equation_of_Time_Gnomical_min = 4 * Equation_of_Time_Gnomical_d

    Longitude_correction_d = Longitude - Time_Zone * 15
    Longitude_correction_min = Longitude_correction_d *4

    EoT_Long_Corrected = (-Equation_of_Time_Gnomical_d + Longitude_correction_d)*4


    # Sun's Altitude & Azimuth

    Observer_True_Hour_Angle_d = (Greanwitch_Mean_Sideral_Time_d + Longitude - Right_Ascension_d)%360
    Observer_True_Hour_Angle_r = radians(Observer_True_Hour_Angle_d)

    Observer_Latitude = radians(Latitude)

    Sun_Altitude_r = asin(sin(Observer_Latitude)*sin(Sun_Declination_r)+cos(Observer_Latitude)*cos(Sun_Declination_r)*cos(Observer_True_Hour_Angle_r))
    Sun_Altitude_d = degrees(Sun_Altitude_r)

    Sun_Zenith_Distance = 90 - Sun_Altitude_d

    local_hour_angle = 15 * ((ST_hrs+(Longitude-Time_Zone*15)/15) - Equation_of_Time_Gnomical_min/60 - 12)
    Sun_Azimuth_r = acos((sin(Sun_Declination_r)-cos(radians(Sun_Zenith_Distance))*sin(Observer_Latitude))/(sin(radians(Sun_Zenith_Distance))*cos(Observer_Latitude)))
    Sun_Azimuth_degree = degrees(Sun_Azimuth_r)

    if (local_hour_angle>0):
        SA_d = 360 - Sun_Azimuth_degree
    else:
        SA_d = Sun_Azimuth_degree

    Sun_Azimuth_d = SA_d % 360
    
    # The Refraction Correction for The Sun's Altitude 

    def atmospheric_pressure():
        P0 = 101325  # Standard atmospheric pressure at sea level (Pa)
        g = 9.81     # Acceleration due to gravity (m/s^2)
        M = 0.0289644  # Molar mass of dry air (kg/mol)
        R = 287.058   # Gas constant for dry air (J/(kg·K))
        T = temperature + 273.15  # Temperature (K)
        return P0 * exp(-g * M * ketinggian_dari_aras_laut / (R * T))

    atmospheric_pressure_milibars= atmospheric_pressure()/100

    def refraction():
        rr = Sun_Altitude_d
        if (rr>15):
            ref=0.00452*tan(radians(Sun_Zenith_Distance))*atmospheric_pressure_milibars/(temperature+273.15)
        else:
            ref = atmospheric_pressure_milibars*(0.1594+0.0196*Sun_Altitude_d + 0.00002 * Sun_Altitude_d**2)/((273.15+temperature)*(1+0.505*Sun_Altitude_d+0.0845*Sun_Altitude_d**2))

        return ref

    humidity= humidityinp/100

    def refractionhumid():
        rr = Sun_Altitude_d
        if (rr>15):
            ref=0.00452*tan(radians(Sun_Zenith_Distance))*atmospheric_pressure_milibars/((1+(7.5e-4*humidity))* (temperature+273.15))
        else:
            ref = atmospheric_pressure_milibars*(0.1594+0.0196*Sun_Altitude_d + 0.00002 * Sun_Altitude_d**2)/((273.15+temperature)*(1+0.505*Sun_Altitude_d+0.0845*Sun_Altitude_d**2)*((1+(7.5e-4*humidity))))
        return ref

    Refraction_d = refraction()
    Refraction_humid = refractionhumid()
    Sun_Altitude_Corrected_d = Sun_Altitude_d - Refraction_humid

    os.system('cls' if os.name=='nt' else 'clear')

    my_banner()
    
    print('Latitude: '+str(Latitude))
    print('Longitude: '+str(Longitude))
    print('Location: ' + location.address)
    print('Date & Time Calculation: ' +Year + '-' + Month + '-' + Second + '' + Hour + ':' + Minute + ':' + Second)
    print('\nGMST: '+str(Greanwitch_Mean_Sideral_Time_d))
    print("Sun's Mean Longitude: " + str(Sun_Mean_Longitude_d))
    print('Perihelion Longitude: ' +str(Perihelion_Longitude_d))
    print('Ecentricity: ' +str(Ecentricity))
    print('Obliquity: ' +str(Obliquity_d))
    print('True Anomaly: ' +str(True_Anomaly))
    print("Sun's True Longitude: " +str(Sun_True_Longitude_d))
    print("Sun's Declination: " +str(Sun_Declination_d))
    print('Right Ascension: ' +str(Right_Ascension_d))
    print('Equation of Time Gnomical in minutes: ' +str(Equation_of_Time_Gnomical_min))
    print('Local Hour Angle: ' +str(local_hour_angle))
    print('Observer True Hour Angle: ' +str(Observer_True_Hour_Angle_d))
    print("Sun's Altitude: " +str(Sun_Altitude_d))
    print("Sun's Zenith Distance: " +str(Sun_Zenith_Distance))
    print("Sun's Azimuth: " +str(Sun_Azimuth_d))
    print('Refraction: ' +str(Refraction_humid))
    print("Sun's Altitude Corrected: " +str(Sun_Altitude_Corrected_d))

    repeat = input("\nDo you want to calculate for other location or other time: (y/n) ")

    if repeat.lower() == "n":
        break
