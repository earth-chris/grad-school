# Assemble URL for download from NASA POWER site

# If you fill in the 40 and -100 for lat and lon, and January 1, 1997 through December 31, 1997 in the date fields at http://power.larc.nasa.gov/cgi-bin/cgiwrap/solar/agro.cgi, then click the "Yes" button next to the question about whether you want ICASA ASCII format, the URL it produces is http://power.larc.nasa.gov/cgi-bin/cgiwrap/solar/agro.cgi?email=agro%40larc.nasa.gov&step=1&lat=40&lon=-100&ms=1&ds=1&ys=1997&me=12&de=31&ye=1997&submit=Yes"

# It should be pretty obvious from staring at the URL that all the parameters (40, -100, etc.) are in there, and what variables they're associated with. Given lat, lon, as well as month, day, and year for both the start and end dates, you can produce any such string. 

# Write a function that takes these parameters and produces the URL string above, just with those values in the correct place. It shouldn't require anything more than the paste() function.

# To download the file, once you have the URL string:
url = "http://power.larc.nasa.gov/cgi-bin/cgiwrap/solar/agro.cgi?email=agro%40larc.nasa.gov&step=1&lat=40&lon=-100&ms=1&ds=1&ys=1997&me=12&de=31&ye=1997&submit=Yes"
txtFile = "myDownloadedFileName.txt"
status = download.file(url, txtFile)
site.df = read.table(txtFile, skip=14)

# Note that it's up to you assign descriptive column names (header=T won't work because it's limited to the case where we just have one row of variable names, followed by tabular data). You can enter your own names this manually, or if you're feeling inspired to learn some more text processing in R, try using readLines() or scan() to get the names and elevation from the preceding lines of the txt file.

# Once you have your function working, you should just be able to write a loop that calls downloadFromURL() once per location, reads into a dataframe, and stores that dataframe somewhere. We recommend a list, though if you're feeling confident with dplyr, rbind'ing all sites into one big dataframe (and including an extra column indicating the site) is also an option. Do not, repeat, not, maintain separate variables for each site's weather dataframe that you would have to operate on individually later.
