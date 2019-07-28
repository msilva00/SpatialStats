library(RCurl)

plotInPlotURL <- getURL("https://raw.githubusercontent.com/luiarthur/rcommon/50ed0069a8d5a33e593796b5c1542f42bca3029c/R/plotInPlot.R", ssl.verifypeer = FALSE)
eval(parse(text = plotInPlot))

myPairsURL = getURL("https://raw.githubusercontent.com/luiarthur/rcommon/50ed0069a8d5a33e593796b5c1542f42bca3029c/R/mypairs.R", ssl.verifypeer = FALSE)
eval(parse(text = myPairsURL))

plotPostsURL = getURL("https://raw.githubusercontent.com/luiarthur/rcommon/50ed0069a8d5a33e593796b5c1542f42bca3029c/R/plotPost.R", ssl.verifypeer = FALSE)
eval(parse(text = plotPostsURL))

colorMultURL = getURL('https://raw.githubusercontent.com/luiarthur/rcommon/50ed0069a8d5a33e593796b5c1542f42bca3029c/R/colorMult.R', ssl.verifypeer = FALSE)
eval(parse(text = colorMultURL))

colorUnderCurveURL = getURL('https://raw.githubusercontent.com/luiarthur/rcommon/50ed0069a8d5a33e593796b5c1542f42bca3029c/R/colorUnderCurve.R', ssl.verifypeer = FALSE)
eval(parse(text = colorUnderCurveURL))