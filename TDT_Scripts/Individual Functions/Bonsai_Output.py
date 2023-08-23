#to take input from the two CSV files and the fMP4 file

import cv2

cap = cv2.VideoCapture("/Users/shivasenthilkumar/Desktop/shiva.mp4")
length_mp4 = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))

with open('/Users/shivasenthilkumar/Desktop/X.csv', 'r') as f1, open('/Users/shivasenthilkumar/Desktop/Y.csv', 'r') as f2:
    fileone = f1.readlines()
    filetwo = f2.readlines()

csv1_len = len (fileone)
csv2_len = len (filetwo)

if (csv1_len == csv2_len  == length_mp4):
    print ("No data was dropped while collecting, all the files (both CSVs and fMPF4 are the same size)")
else:
    print ("Oh no! Some Data was dropped, look through all the file lengths again!")
    print ("First CSV" + csv1_len)
    print ("Second CSV" + csv2_len)
    print ("FMP4" + length_mp4)
