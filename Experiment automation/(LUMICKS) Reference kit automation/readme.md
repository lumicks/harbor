<h1>Automation script for reference kit</h1><br>
<b>Author</b>: Maartje Pontier, Aafke van den Berg<br>
<b>Key words</b>: kymograph, fluorescence, reference kit<br>
<b>Research subjects</b>: DNA-binding proteins<br>
<br>
<p>This script automates the workflow for the LacI reference kit. 
<p>The script follows the following steps: <br></p><div>1) Catch beads<br>2) Catch a DNA tether<br>3) Check if it is only a single tether, if not, try to break tethers until there is 1 left.<br>4) Move microstage to protein channel <br> 5) Stretch the tether to a force of 10 pN<br>6) Start a kymograph <br>7) Stop the kymograph <br> 9) Repeat steps 1 to 7 </div><br>

The script has been tested rigourously, but may need be optimized for your system.
Parameters that you, for example, may have to optimize are: the pressure of the microfluidics for obtaining optimal flow, the names of the channels in the flowcell and the number of kymographs that you want to record.
<br>
</em></p><p><br></p><div class="se-component se-image-container __se__float-none"><figure style="margin: 0px;"><img style="" data-origin="," data-file-size="31661" data-file-name="kymo_10_B23.png" data-percentage="auto,auto" data-align="none" data-size="," data-rotatey="" data-rotatex="" data-proportion="true" data-rotate="" alt="" src="kymo_10_B23.png" data-index="3"></figure></div><p><em>Figure 1: Kymograph recorded using the LacI reference kit.</em></p><p><br></p><p><em>
<br>
