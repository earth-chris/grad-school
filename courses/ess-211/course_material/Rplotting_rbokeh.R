# http://hafen.github.io/rbokeh/

library(rbokeh)
library(plyr) # will be covered in R session 6


# --- sample data ---
data(airquality) # loads the airquality dataset into the workspace
attach(airquality) # extracts each variable as a vector

# Add month abbreviation
airquality$Mon <- factor(month.abb[airquality$Month], levels=month.abb[unique(airquality$Month)])


# ----- Simple scatter plot ----- 
p <- figure(title="Ozone vs. Temperature in NY, 1973") %>%   # Initialize a Bokeh figure
  ly_points(Temp, Ozone, data=airquality,                    # Add a “points” layer to a Bokeh figure 
            glyph = 17,                                      # glyph type for points
            hover = list(Temp, Ozone))  # when a point is hovered over in the plot, a list of the specified variable values will be shown 

p                                                            # Display the figure


# ----- Simple line plot -----  
p <- figure() %>%                           # Initialize a Bokeh figure
  ly_lines(Temp, data=airquality)           # Add a “lines” layer to a Bokeh figure 

p                                           # Display the figure


# ----- Simple point and line plot ----- 
p <- figure() %>%                                       # Initialize a Bokeh figure
  ly_lines(Ozone, data=airquality) %>%                  # Add a “lines” layer to a Bokeh figure
  ly_points(Temp, data=airquality, color='blue') %>%    # Add a “points” layer to a Bokeh figure
  ly_lines(Temp,  data=airquality, color='blue') %>%    # Add a “lines” layer to a Bokeh figure 
  ly_lines(Solar.R,  data=airquality, color='red')      # Add a “lines” layer to a Bokeh figure

p                                                       # Display the figure


# ----- Histogram plot -----
h <- figure() %>%                                       # Initialize a Bokeh figure
  ly_hist(Temp, data=faithful, breaks=50)               # Add a “histogram” layer to a Bokeh figure

h                                                       # Display the figure


# ----- Box plot -----
b <- figure(title="Boxplot of Ozone values") %>%        # Initialize a Bokeh figure
  ly_boxplot(y=Ozone, x=Mon, data=airquality)           # Add a “boxplot” layer to a Bokeh figure

b                                                       # Display the figure


# ----- Bar plot -----
# Calculate monthly average Temperature
monthly_averages <- ddply(airquality, .(Mon), summarise, avg = mean(Temp), desc='Temp')

b <- figure(title="Monthly Average Temperature") %>%          # Initialize a Bokeh figure
  ly_bar(Mon, avg, data=monthly_averages) %>%                 # Add a “bar plot” layer to a Bokeh figure
  ly_text(Mon, avg, text=round(avg,1), data=monthly_averages) # Add a "text" layer to a Bokeh figure

b                                                             # Display the figure


# ------ Two bar plots dodged -------
# Calculate monthly average Ozone
# ddply function will be covered in R session 6
tmp              <- ddply(airquality, .(Mon), summarise, avg=mean(Ozone,na.rm=T), desc='Ozone')
# Combine the average temperature and ozone data
monthly_averages <- rbind(monthly_averages, tmp)

b <- figure(title="Monthly Average Temperature") %>%                    # Initialize a Bokeh figure
  ly_bar(Mon, avg, color=desc, data=monthly_averages, position='dodge') # Add a “bar plot” layer to a Bokeh figure

b                                                                       # Display the figure


# ------ Two bar plots stacked -------
b <- figure(legend_location = "top_left") %>%                           # Initialize a Bokeh figure
  ly_bar(desc, avg, color=Mon, data=monthly_averages, position='stack') # Add a “bar plot” layer to a Bokeh figure

b                                                                       # Display the figure
