# R Code and data for the MS: "Behavioural thermoregulation and microclimate reshape climate-driven range forecasts"
Juan G. Rubalcaba*, Guillermo Fandos, José A. Díaz
Department of Biodiversity, Ecology and Evolution. Faculty of Biological Sciences. Complutense University of Madrid
*jg.rubalcaba@gmail.com

Code and data supporting analyses of behavioural buffering and hybrid species distribution models under climate warming.
### Abstract
Microclimate heterogeneity and behavioural thermoregulation can buffer ectotherms against environmental warming, yet these processes are rarely captured in models forecasting species’ responses to climate change. Here we use microclimate and biophysical models validated against field data separated by 25 years of global warming (1997-2022) to quantify how behavioural thermoregulation modifies body temperature in a widespread lizard species, and how this influences its geographic distribution. We show that body temperatures increased by only 18-29% of the rise in environmental temperatures, showing the role of behavioural buffering. When incorporated into species distribution models, this buffering mechanism emerged as a key predictor of current distributions and led to markedly different projections under future warming compared to conventional climate-based models. These differences revealed cryptic refugia in cooler regions, where behavioural buffering reduces climate impacts, and false refugia in warmer regions, where buffering capacity is exceeded and persistence is overestimated by correlative approaches. Our results show that accounting for how organisms interact with microclimates can fundamentally alter forecasts of species’ responses to climate change, highlighting the importance of integrating behavioural and mechanistic processes into predictive models.

## Repository content: 
<p> <b> RCode_biophysical_model.R: </b>  Code to generate mechanistic layers (thermoregulatory inaccuracy and thermoregulatory window) </p> 
<p> <b> Sources file </b> </p>
    <p> <b> SIMULATION RESULTS </b> </p>
    <p> - dataMAY.R: Body temperatures, operative temperatures, inaccuracy and thermoregulation window in May (current) </p>
    <p> - dataJUNE.R: Body temperatures, operative temperatures, inaccuracy and thermoregulation window in June (current) </p>
    <p> - dataMAY_warm.R: Body temperatures, operative temperatures, inaccuracy and thermoregulation window in May (future) </p>
    <p> - dataJUNE_warm.R: Body temperatures, operative temperatures, inaccuracy and thermoregulation window in June (future) </p>
    <p> <b> RASTER LAYERS </b>  </p>
    <p> - map.grd/.gri: base raster layer </p>
    <p> - buffer_map.grd/.gri: predicted buffer (delta Tb/delta Te) </p>
    <p> - meandb_map.grd/.gri: thermoregulatory inaccuracy (current conditions) </p>
    <p> - meanActivity_map.grd/.gri: thermoregulatory window (current conditions) </p>
    <p> - meandb_warm_map.grd/.gri: thermoregulatory inaccuracy (future) </p>
    <p> - meanActivity_warm_map.grd/.gri: thermoregulatory window (future) </p>
    
# Session Info
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26200)

Attached packages:
 [1] ggplot2_3.5.2      RColorBrewer_1.1-3 RNetCDF_2.9-2      RNCEP_1.0.10       maps_3.4.2.1       terra_1.8-60      
 [7] lubridate_1.9.3    stringr_1.5.1      dplyr_1.2.0        tidyr_1.3.1        raster_3.6-32      sp_2.2-0          
[13] microclima_0.1.0   NicheMapR_3.3.2  
