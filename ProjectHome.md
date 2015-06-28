# JS Planar Ordinary Voronoi Tesselation #

created by Greg Ross

---


This is a pure JavaScript implementation of an incremental algorithm for a planar ordinary Voronoi diagram. The algorithm uses a quaternary tree for spatial indexing to enhance performance and facilitate efficient nearest neighbor searches. The geometry is represented by a winged-edge data structure. This provides a succinct representation from which information such as the minimum spanning tree can be derived.

An incremental algorithm was chosen over Fortune's sweepline solution so that new points can be inserted into the diagram without needing to regenerate the whole tessellation. The use of the quaternary tree also reduces the average time complexity to **O(n)**.

Tested in Chrome, Firefox, Opera, Safari and Internet Explorer 6. **Works in IE via use of excanvas.** (See test.html in the download.)

## See [here](http://www.grvisualisation.50webs.com/javascript_voronoi.html) for a working example. Also, see the [jsFiddle](http://jsfiddle.net/gregross/dcdw2/) example. ##

## See [here](http://javascript-voronoi.googlecode.com/svn/trunk/googleVizApi.html) for API details. ##

![http://javascript-voronoi.googlecode.com/svn/trunk/images/voronoiScreenShot.png](http://javascript-voronoi.googlecode.com/svn/trunk/images/voronoiScreenShot.png)

Since VML renders very slowly in IE it is advised that the `fillPolygons` option is set to **false** when the script is run in IE, especially when there are a substantial number of points in the diagram.

## Features ##

  * **fast** incremental algorithm: potentially faster than Fortune's
  * add points without entire redraw
  * **pure JavaScript** implementation. **No need for Flash**.
  * supports **tooltips** over Voronoi regions
  * customisable **colour gradients**
  * works in all popular browsers
  * wrapped in the **Google Visualization API**



## Applications ##

  * spatial distributions, e.g. animal territories and bird nesting patterns
  * nearest neighbor search
  * terrain modelling
  * robot path planning
  * Delaunay Triangulation and mesh generation

<br>

The following book provides an excellent description of the above algorithm. This is what was followed in its implementation:<br>
<br>
Okabe A., Boots B., Sugihara K., Chiu S.N.: "Spatial Tessellations. Concepts and Applications of Voronoi Diagrams." Wiley, Chichester, second edition (2000).<br>
<br>
<br>

<hr />

Other projects by Greg Ross...<br>
<br>
<br>
<h3><i><b><a href='http://www.grvisualisation.50webs.com/visi_scroll.html'>Visi-Scroll</a>.</b></i></h3>

<b>An HTML scrollbar with highligting for easy access to information</b>


<br>

<h3><i><b><a href='http://www.grvisualisation.50webs.com/excelrangefinder.html'>Excel RangeFinder</a>.</b></i></h3>

<b>An overview+detail visualisation for navigation in Excel spreadsheets</b>


<br>

<h3><i><b><a href='http://www.grvisualisation.50webs.com/'>Magic Table</a>.</b></i></h3>

<b>A JavaScript library that allows you to see more in your data by applying some simple visual techniques to transform a table.</b>


<br>

<h3><i><b><a href='http://code.google.com/p/webgl-surface-plot/'>WebGL surface plot</a>.</b></i></h3>

<b>3D surface plotting in JavaScript via WebGL</b>


<br>

<h3><i><b><a href='http://code.google.com/p/javascript-surface-plot/'>JavaScript 3D surface plotting</a>.</b></i></h3>

<b>A pure JavaScript implementation of a 3D surface plot</b>


<hr />

<a href='http://itunes.apple.com/gb/app/claptrax/id427145886?mt=8'><img src='http://www.grvisualisation.50webs.com/images/clapForMusic/SplashViewPhone.png' /></a>

<hr />

<a href='http://ax.itunes.apple.com/gb/app/peekaboo-3d/id380735126?mt=8'><img src='http://www.grvisualisation.50webs.com/images/peekaboo/PeekabooAd.png' /></a>