<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>

    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>




<!-- GETTING STARTED -->
<h2 id="getting-started">Getting Started</h2>
<p>This is my homework for the linear algebra part of TC course. This README file should help you set up the code, build all executables and execute them in one simple step.</p>
<h3 id="prerequisites">Prerequisites</h3>
<p>To run the script &quot;compile.sh&quot; you will need installed bash on your system and CMake program for building software from source files and managing dependencies. To compile the source code, GNU Fortran compiler is recommended, but other compilers should also work. In case the &quot;compile_and_run.sh script isn't working, you can always manually compile source files. To plot all the scripts you will need Pytyhon 3.5+ with numpy, matplotlib and energydiagram modules installed, and GNU Octave and GNU Plot programs. </p>
<ul>
<li><p>GNU Fortran</p>
<pre><code>
</code></pre>
</li>
<li><p>CMake Version 3.5+ (optional)</p>
<pre><code>
</code></pre>
</li>
<li><p>Python 3.5+ with numpy, matplotlib and energydiagram modules (optional)</p>
<pre><code>
</code></pre>
</li>
<li><p>GNU Octave and GNU Plot (optional)</p>
<pre><code>
</code></pre>
</li>
</ul>
<p>Below is an example how to install everything on Linux Ubunutu. If you are using other linux distribution or MacOS, you can find the information online.</p>
<ol>
<li><pre><code class="lang-sh">sudo apt <span class="hljs-keyword">install</span> gfortran
</code></pre>
</li>
<li><pre><code class="lang-sh">sudo apt <span class="hljs-keyword">install</span> cmake
</code></pre>
</li>
<li><pre><code class="lang-sh">sudo apt <span class="hljs-keyword">install</span> gnuplot
</code></pre>
</li>
<li><pre><code class="lang-sh">sudo apt <span class="hljs-keyword">install</span> octave
</code></pre>
</li>
<li><pre><code class="lang-sh">pip <span class="hljs-keyword">install</span> numpy, matplotlib, energydiagram
</code></pre>
</li>
<li>If everything went okay (with no errors) you should now have directory bin with all executables.<p align="right">(<a href="#readme-top">back to top</a>)</p>



</li>
</ol>
<!-- USAGE EXAMPLES -->
<h2 id="usage">Usage</h2>
<p>Source files are consisted of 4 modules:<br>linear_algebra - This module constains all functions and procedures related for the math manipulations of the Part 1 of the task.<br> math - This module contains necessary procedures for math manipulations for the Part 2 of the task.<br> kinds - This module defines kinds of real numbers used to ensure the same precision in every machine and compiler. <br> inout - This module contains procedures for the input and the output of the files.<br>The functions and procedures are generally well-documented.  <br/><br/> If you have any question, you can send me an email.</p>
<p>After Running ./compile_and_run.sh you should have all the output files nicely presented in different folders. In every folder with output files you will find scripts that create plots that you can see in the pdf file. If you have all the required packages instaleld and you agreed to plot while running ./compile_and_run.sh you should also have the plots there.</p>
<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ROADMAP -->
<h2 id="roadmap">Roadmap</h2>
<ul>
<li>bin - This is where are all the binaries </li>
<li>src - This is where are all the source files</li>
<li>INPUT - Folder where are the input files.</li>
<li>OUTPUT - Folder where are the output files.</li>
<li>SCRIPTS - Folder where are the python/octave scripts for plotting.</li>
<li>CMakeLists.txt - Don&#39;t touch this file. It configures installation with CMAKE</li>
<li>README<p align="right">(<a href="#readme-top">back to top</a>)</p>




</li>
</ul>
<!-- CONTACT -->
<h2 id="contact">Contact</h2>
<p>Zarko Ivkovic  zivkoviv7@alumnes.ub.edu</p>
<p>Please email me if you have troubles configuring the installation or if you have any questions about code that you want to discuss.</p>
<p align="right">(<a href="#readme-top">back to top</a>)</p>

