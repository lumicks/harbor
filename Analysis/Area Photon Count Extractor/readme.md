<h1>Area Photon Count Extractor</h1><br>
<b>Author</b>: John Watters<br>
<b>Key words</b>: fluorescence<br>
<b>Research subjects</b>: DNA-binding proteins, Other applications<br>
<br>
<p>This script is used to extract the sum of photon counts in a line scan from chosen areas of a kymograph in a .h5 file. This analysis method is useful when you want to measure raw intensity values of fluorescent molecules to quantify how many molecules are present without worrying about potential biases introduced by using line tracking methods.</p><p>You can download the script <a href="https://github.com/watters16j/CTrap-Scripts">here</a>.</p><p>&nbsp;</p><p>The user will be asked via text inputs to navigate to the .h5 file of interest and apply a multiplier value to photon counts if desired (for display purposes only monitoring the method used in CTrapVis.py). The user will then choose how they want to define the area to analyze either by [1] Explicitly defined dimensions or [2] dragging a window on the kymograph in a pop-up window. Then the user will be asked to click on other regions to extract where the click defines the center-left point of the box (in the cartoon below,&nbsp; '-' is the number of time points being analyzed, '|' is the number of pixels being analyzed, and 'X' marks the spot where the user clicks).</p><p>&nbsp;</p><p>--------------</p><p>|&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;|</p><p>|&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;|</p><p>X&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; |</p><p>|&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;|</p><p>|&nbsp; &nbsp; &nbsp;&nbsp; &nbsp; &nbsp; &nbsp; &nbsp;|</p><p>--------------</p><p>&nbsp;</p><p>Once the points are selected on the plot, then you can exit the plot. This will show you a separate image showing the regions the script will extract. If you used the&nbsp;'dragging a kymograph window' method of defining the area dimensions, then you will see the original box drawn in blue with the additionally clicked boxes drawn in orange.</p><p><br></p><div class="se-component se-image-container __se__float-none"><figure style="margin: 0px;"><img src="img0.png" alt="" data-rotate="" data-proportion="true" data-rotatex="" data-rotatey="" data-size="," data-align="none" data-percentage="auto,auto" style="" data-file-name="showing_regions_to_extract.png" data-file-size="958215" data-origin="," data-index="0"></figure></div><p><br></p><p>If the user is happy with the point selection, they can confirm that through user input and the script will be extracted the data to a .xlsx file. There are two types of sheets in the output .xlsx. The first sheet-type shows the sum of the columns of each region (at each time point) and the extracted simple statistics from each region (average and standard deviation) for each channel you want to extract.The second sheet-type records some metadata that might be useful if you need to remake any plots/redo any analysis. In addition, an image of the extracted regions is also outputted as a .png for a quick resource on which region #'s align with other regions of the plot.</p><p><u>Example .png output</u></p><p><br></p><div class="se-component se-image-container __se__float-none"><figure style="margin: 0px;"><img src="img1.png" alt="" data-rotate="" data-proportion="true" data-rotatex="" data-rotatey="" data-size="," data-align="none" data-percentage="auto,auto" data-file-name="20220306-150211 Kymograph 16_extracted_intensity_regions_06_12_2022.png" data-file-size="187417" data-origin="," style="" data-index="1"></figure></div><p></p><p><u>Example .xlsx output</u></p><p>Example Photon Count Data</p><div class="se-component se-image-container __se__float-none"><figure style="margin: 0px;"><img src="img2.png" alt="" data-proportion="true" data-align="none" data-file-name="DzhPDCfPRnB4AAAAAElFTkSuQmCC" data-file-size="0" data-origin="," data-size="," data-rotate="" data-rotatex="" data-rotatey="" data-percentage="auto,auto" style="" data-index="2"></figure></div><p>Example Metadata</p><div class="se-component se-image-container __se__float-none"><figure style="margin: 0px;"><img src="img3.png" alt="" data-proportion="true" data-align="none" data-file-name="ATjPX3cdhevKAAAAAElFTkSuQmCC" data-file-size="0" data-origin="," data-size="," data-rotate="" data-rotatex="" data-rotatey="" data-percentage="auto,auto" style="" data-index="3"></figure></div><p><br></p><p><strong><u>Version Information:</u></strong></p><ul type="disc"> <li>Lumicks.pylake - 0.8.1</li><li>pandas - 1.4.1</li><li>matplotlib - 3.5.1</li><li>numpy - 1.22.3</li></ul><p>Enjoy!</p><p>&nbsp;</p><br><br>
<b>If you like this script, please cite me</b>: <p>Watters, J.W. (2022) Area Photon Count Extractor. Retrieved from https://harbor.lumicks.com/<br></p><br>
<br>
