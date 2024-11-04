<h1>Catching and generating ds-ssDNA hybrids</h1><br>
<b>Author</b>: Olivia Yang<br>
<b>Key words</b>: automated script, WLC, dna, Overstretching, automation, DNA, script<br>
<b>Research subjects</b>: DNA-binding proteins<br>
<br>
<p>This script automates the process of</p><ol><li>catching beads</li><li>catching DNA - user specifies length (nucleotides) @ line 38, variable total_nt<br></li><li>overstretching to generate ssDNA - user specifies how many nucleotides to unwind @ line 39, variable ssDNA_nt. If no unwinding required, set as 0.<br></li><li>checks distance against expected model eWLC+FJC for hybrid DNA (if unwound).</li><li>ends in buffer channel with flow off and at low tension<br></li></ol><p>This script expects the user has<br></p><ul><li>caught beads, set template, and calibrated force</li><li>beads in channel 1, DNA in channel 2, low salt buffer in channel 3</li><li>DNA to generate hybrid ss-dsDNA has nicks to unwind specified amount<br></li></ul><p>Contributors:</p><ul><li>Roeland van Wijk</li><li>Olivia Yang<br></li></ul><br><br>
<br>

