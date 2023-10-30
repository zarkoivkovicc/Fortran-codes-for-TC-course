<!-- GETTING STARTED -->
<h2 id="getting-started">Getting Started</h2>
<p>This is my homework for the calculus part of TC course. This README file should help you set up the code and build all executables in one simple step.</p>
<h3 id="prerequisites">Prerequisites</h3>
<p>To run the script &quot;compile.sh&quot; you will need installed bash on your system and CMake program for building software from source files and managing dependencies. To compile the source code, GNU Fortran compiler is recommended, but other compilers should also work. In case the compile.sh script isn't working, you can always manually compile source files. </p>
<ul>
<li><p>bash -
Bash usually comes pre-installed in all linux distros</p>
</li>
<li><p>CMake Version 3.5+</p>
<pre><code>
</code></pre><h3 id="installation">Installation</h3>
</li>
</ul>
<p>Below is an example how to set up everything.</p>
<ol>
<li>If you don't have CMake installed, you can install it on Ubuntu using the command:<pre><code class="lang-sh">sudo apt <span class="hljs-keyword">install</span> cmake
</code></pre>
</li>
<li>If you are using other linux distributions, the information how to install cmake can be found online<pre><code class="lang-sh">
</code></pre>
</li>
<li>Run this command <pre><code class="lang-js">./compile.sh
</code></pre>
</li>
<li>If everything went okay (with no errors) you should now have directory bin with all executables.<p align="right">(<a href="#readme-top">back to top</a>)</p>



</li>
</ol>
<!-- USAGE EXAMPLES -->
<h2 id="usage">Usage</h2>
<p>Source files are consisted of 3 modules:<br>integrals - This module constains all functions and procedures related for the Part 1 of the task.<br> oprimization - This module contains necessary    procedures for the Part 2 of the task.<br> kinds - This module defines kinds of real numbers used to ensure the same precision in every machine and compiler.<br>The functions and procedures are generally well-documented.  <br/><br/> If you have any question, you can send me an email.</p>
<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ROADMAP -->
<h2 id="roadmap">Roadmap</h2>
<ul>
<li>bin - This is where all the binaries are  </li>
<li>src - This is where all source files are</li>
<li>data - Folder where are the examples of output files(if they exist).</li>
<li>CMakeLists.txt - Don&#39;t touch this file. It configures installation with CMAKE</li>
<li>README<p align="right">(<a href="#readme-top">back to top</a>)</p>




</li>
</ul>
<!-- CONTACT -->
<h2 id="contact">Contact</h2>
<p>Zarko Ivkovic  zivkoviv7@alumnes.ub.edu</p>
<p>Please email me if you have troubles configuring the installation or if you have any questions about code that you want to discuss.</p>
<p align="right">(<a href="#readme-top">back to top</a>)</p>
