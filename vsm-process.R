# Options -----------------------------------------------------------------

display.mh = "no"
export.data = "yes"
export.plots = "yes"
export.grid = "yes"


# Import packages ---------------------------------------------------------

my.packages <- c("rChoiceDialogs", "purrr", "dplyr", "ggplot2", "ggpubr", "scales", "sfsmisc", "signal")
lapply(my.packages, require, character.only = TRUE)


# Functions ---------------------------------------------------------------

shift.data <- function(data) {
  moment.max = max(data)
  moment.min = min(data)
  moment.diff = abs(moment.max) - abs(moment.min)
  if ((moment.max < 0) |
      (moment.min > 0)) {
    stop("Data shift not possible, bad data")
  }
  shifted.data = data - (moment.diff / 2)
  return(shifted.data)
}


read.dat <- function(flnm) {
  read.csv(flnm, skip = 12, sep = "\t", header = F, col.names = c("field", "moment", "blank")) %>%
    mutate(filename = flnm) %>%
    select(field, moment, filename) %>%
    na.omit() %>%
    group_by(filename) %>%
    mutate(range = round((max(field)-min(field))/2,0)) %>%
    ungroup() %>%
    select(range, field, moment)
}

read.conc <- function(flnm) { read.csv(flnm, header = FALSE, skip = 0) }


find.fit <- function (data, type){
  tolerance <- 0.9999
  tolerance.drop <- 0.0001
  
  if (type == "middle"){
    frac <- 0.2
    data <- data %>% 
      dplyr::filter(between(moment, frac * min(moment), frac * max(moment))) %>% 
      arrange(field)
    # plot(data$field,data$moment)
  }
  
  if (type == "negative.reciprocal") {
    frac <- 0.8
    data <- data %>% 
      dplyr::filter(between(moment, min(moment), frac * min(moment))) %>% 
      arrange(field)
  }
  
  if (type == "positive.reciprocal") {
    frac <- 0.8
    data <- data %>% 
      dplyr::filter(between(moment, frac * max(moment), max(moment))) %>% 
      arrange(field)
  }
  
  fits = data.frame()
  
  while (tolerance > 0.75) {
    i = k = 1
    j = nrow(data)
    
    while (i < j & j - i > 8) {
      if (nrow(fits) < k){
        if (type == "middle") {
          fit = lm(data$moment[i:j] ~ data$field[i:j])
          summary(fit)$coefficients[2]
          fits[k, 1] = type
          fits[k, 2] = data$field[i]
          fits[k, 3] = data$field[j]
        }
        
        if (type == "negative.reciprocal") {
          fit = lm(data$moment[i:j] ~ data$reciprocal[i:j])
          fits[k, 1] = type
          fits[k, 2] = data$reciprocal[i]
          fits[k, 3] = data$reciprocal[j]
        }
        
        if (type == "positive.reciprocal") {
          fit = lm(data$moment[i:j] ~ data$reciprocal[i:j])
          fits[k, 1] = type
          fits[k, 2] = data$reciprocal[i]
          fits[k, 3] = data$reciprocal[j]
        }
        
        fits[k, 4] = j - i + 1
        fits[k, 5] = summary(fit)$coefficients[2]
        fits[k, 6] = summary(fit)$coefficients[1]
        fits[k, 7] = summary(fit)$r.squared
      } 

      if (fits[k, 7] > tolerance) {
        colnames(fits) = c("source" ,"lower.bound", "upper.bound", "num.points", "slope","intercept", "r.2")
        return(fits[k,])
      }
    
      if (type == "middle") {    
        i = i + 1
        j = j - 1
      }
      
      if (type == "negative.reciprocal") {    
        j = j - 1
      }
      
      if (type == "positive.reciprocal") {    
        i = i + 1
      }
      k = k + 1
    }

    tolerance = tolerance - tolerance.drop
    # print(tolerance)
  }
  stop("Failed to meet fit tolerance")
}

# Function to calculate the nanoparticle size based on Chantrell fitting
calc.d <- function(kB, temperature, Xi, Ho, moment, magnetization) {
  diameter = ((((18 * ( kB ) * temperature) / (pi * magnetization)) * sqrt(Xi / (3 * moment * Ho))) ^ (1 / 3)) * (10 ^ 9)
  return(diameter)
}


# Function to calculate the nanoparticle size distribution based on Chantrell fitting
calc.sigma <- function(Xi, Ho, moment) {
  sigma = ((log(3 * Xi / (moment / Ho))) ^ (1 / 2)) / 3
  return(sigma)
}


# Function to calculate 1st standard deviation
calc.std.dev <- function(diameter, sigma) {
  std.dev.low = diameter / exp(sigma)
  std.dev.high = diameter * exp(sigma)
  std.dev = c(std.dev.low, std.dev.high)
  return(std.dev)
}


# Function to calculate the real saturation magnetization of the sample
calc.Ms <- function(moment, concentration.fe, volume, density) {
  concentration.fe3o4 = concentration.fe / 0.72
  Am2 = moment
  Am2mg = Am2 / (concentration.fe3o4 * volume / 1000)
  Am = Am2mg * density
  kAm = Am / 1000
  return(kAm)
}

# Function to plot log normal distribution on exportable plot
log.normal = function(x, size, sigma) {
  (1 / ((x / size) * sigma * sqrt(2 * pi))) * exp((-log(x / size) ^ 2) / (2 * sigma ^ 2))
}


# Import data -------------------------------------------------------------

setwd("/Users/ericteeman/Google Drive/Research/Data/VSM/")
setwd(rchoose.dir(caption = "Select Directory")) # Asks user to choose directory containing data files

dat <- list.files(pattern = "*\\d.txt", full.names = T, recursive = F) %>% map_df(~ read.dat(.))


# Spectrometer information ------------------------------------------------

# Set physical values and constants
vol <- 100 #uL
density <- 5180 # kg/m^3
kB <- 1.38e-23 # J/K
temperature <- 298 # K


if (file.exists("conc.txt")) {
  conc <- list.files(pattern = "conc.txt", full.names = T, recursive = F) %>% 
    map_df(~ read.conc(.))
  conc <- conc[, 1] #mgFe/mL
  conc.new = conc / 1000 #gFe/mL
  vol.new = vol / 1000 #mL
  mass = conc.new * vol.new #gFe
  mh.label = expression(paste("M [kA m" ^ "-1", "]"))
} else {
  conc <- 0 #mgFe/mL
  mass <- 1 #filler value to prevent conversion without known conc
  mh.label = expression(paste("m [Am" ^ "2", "]"))
}

# Data correction to account for instrument error
dat = dat %>%
  group_by(range) %>%
  mutate(moment = shift.data(moment)) %>% #center moment around zero
  ungroup() %>%
  mutate(range = range * (10 ^ -1)) %>% #convert gauss to mT
  mutate(field = field * (10 ^ -1)) %>% #convert gauss to mT
  mutate(reciprocal = 1 / field) %>%
  mutate(moment = moment * (10 ^ 3)) %>% #convert emu to Am2
  mutate(norm = 2 * ((moment - min(moment)) / (max(moment) - min(moment))) - 1) %>%
  mutate(magnetization = ((moment / ((conc / 0.72) * vol / 1000)) * density) / 1000) #calculate magnetization

fits = find.fit(dat, "middle") %>%
  dplyr::bind_rows(., find.fit(dat, "negative.reciprocal")) %>%
  dplyr::bind_rows(., find.fit(dat, "positive.reciprocal")) %>%
  group_by(source) %>%
  mutate(Xi = abs(slope)) %>%
  mutate(Ho = abs(1 / (-intercept/slope))) %>%
  mutate(moment = abs(intercept)) %>%
  mutate(Ms = case_when(conc == 0 ~ 446,
                        conc != 0 ~ abs(calc.Ms(moment, conc, vol, density))))

# View(fits)

info = data.frame(
  conc,
  round(mean(c(fits$Ms[2], fits$Ms[3])), 2),
  round(calc.d(kB, temperature, fits$Xi[1],
         mean(c(fits$Ho[2],fits$Ho[3])),
         mean(c(fits$moment[2], fits$moment[3])),
         mean(c(fits$Ms[2], fits$Ms[3]))), 2),
  round(calc.sigma(fits$Xi[1], mean(fits$Ho[2],fits$Ho[3]), mean(fits$moment[2],fits$moment[3])), 2))
colnames(info) = c("Conc [mgFe/mL]", "Ms [kA/m]", "Size [nm]", "Sigma")


# Plots -------------------------------------------------------------------

theme_new <- function (base_size=24, base_line_size=1) {
  theme_bw(base_size=base_size,
           base_family="") %+replace%
    theme(
      axis.text.x = element_text(size=base_size, margin = margin(t=0.75*base_size,b=0.25*base_size), color="black"),
      axis.text.y = element_text(size=base_size, margin = margin(r=0.75*base_size,l=0.25*base_size), color="black"),
      axis.title=element_text(size=base_size, color="black"),
      axis.line=element_line(size=base_line_size, lineend="square", color="black"),
      axis.ticks=element_line(size=base_line_size, lineend="square", color="black"),
      axis.ticks.length=unit(-8, "pt"),
      panel.border=element_blank(),
      panel.grid=element_blank(),
      aspect.ratio = 1,
      legend.background=element_rect(fill="transparent", colour=NA),
      legend.position=c(1, 1),
      legend.justification=c("right", "top"),
      legend.direction="vertical",
      legend.title=element_text(size=0.75*base_size),
      legend.title.align=0.5,
      legend.text=element_text(size=0.75*base_size),
      legend.text.align=0
    )
}

xlab = expression(paste(mu[0],"H [mT]"))
ylab = mh.label
# ylab = "Norm. [a.u.]"

data.set <- dat %>% dplyr::filter(range == min(range))

p1 = ggplot(data.set, aes(x = field, y = magnetization)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_path(size = 1) +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  scale_y_continuous(breaks = c(min(pretty(min(data.set$magnetization))), 0, max(pretty(max(data.set$magnetization)))), 
                     limits = c(min(pretty(min(data.set$magnetization))), max(pretty(max(data.set$magnetization))))) +
  theme_new() +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)

data.set <- dat %>% dplyr::filter(range == median(range))

p2 = ggplot(data.set, aes(x = field, y = magnetization)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_path(size = 1) +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  scale_y_continuous(breaks = c(min(pretty(min(data.set$magnetization))), 0, max(pretty(max(data.set$magnetization)))), 
                     limits = c(min(pretty(min(data.set$magnetization))), max(pretty(max(data.set$magnetization))))) +
  theme_new() +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)

data.set <- dat %>% dplyr::filter(range == max(range))

p3 = ggplot(data.set, aes(x = field, y = magnetization)) +
  geom_point(size = 3, shape = 1, stroke = 1) +
  geom_path(size = 1) +
  scale_x_continuous(breaks = pretty_breaks(3)) +
  scale_y_continuous(breaks = c(min(pretty(min(data.set$magnetization))), 0, max(pretty(max(data.set$magnetization)))), 
                     limits = c(min(pretty(min(data.set$magnetization))), max(pretty(max(data.set$magnetization))))) +
  theme_new() +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)


data.set = data.frame(x = 0)
xmin = 0.5*info$`Size [nm]`
xmax = 1.5*info$`Size [nm]`
xlab = expression(paste(d[0]," [nm]"))
ylab = "Density [%]"

p4 = ggplot(data.set) +
  stat_function(fun = log.normal, args = list(info$`Size [nm]`, info$Sigma), geom="line", size = 1) +
  scale_x_continuous(breaks = pretty_breaks(3), limits = c(xmin, xmax)) +
  scale_y_continuous(breaks = pretty_breaks(3)) +
  theme_new() +
  guides(col = guide_legend(ncol = 1)) +
  labs(x = xlab, y = ylab)


# Set export directory whether or not saving images is selected

if (export.data == "yes" || export.plots == "yes" || export.grid == "yes") {
  main.directory = getwd()
  export.directory = "export"
  dir.create(file.path(main.directory, export.directory), showWarnings = FALSE)
  setwd(file.path(main.directory, export.directory))
}

sc = 1    # Set scaling for saved images

if(export.data == "yes"){
  write.csv(info, "info.csv", row.names = FALSE)
  write.csv(dat, "dat.csv", row.names = FALSE)
}

if(export.plots == "yes"){
  ggsave("mh.png", p2, width = 4.5 * sc, height = 4.5 * sc, dpi = "retina")
  ggsave("histogram.png", p4, width = 4.5 * sc, height = 4.5 * sc, dpi = "retina")
}

if(export.grid == "yes"){
  grid <- ggarrange(p1, p2, p3, p4, hjust = -0.25, heights = c(1,1), widths = c(1,1),
                        labels = c("(a)", "(b)","(c)","(d)"), font.label = list(size = 24),
                        ncol = 2, nrow = 2)
  ggsave("grid.png", grid, width = 9 * sc, height = 9 * sc, dpi = "retina")
}

if (display.mh == "yes"){
  p1
}

info
