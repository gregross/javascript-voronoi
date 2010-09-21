/*
 * Register the name space
 * ***********************
 */
function registerNameSpace(ns)
{
    var nsParts = ns.split(".");
    var root = window;
	var n = nsParts.length;
    
    for (var i = 0; i < n; i++) 
    {
        if (typeof root[nsParts[i]] == "undefined") 
            root[nsParts[i]] = new Object();
        
        root = root[nsParts[i]];
    }
}

registerNameSpace("greg.ross.visualisation");

/*
 * This is the main class and entry point of the tool
 * and represents the Google viz API.
 * ***************************************************
 */
greg.ross.visualisation.Voronoi = function(container)
{
    this.containerElement = container;
}

greg.ross.visualisation.Voronoi.prototype.draw = function(data, options)
{
	var xPos = options.xPos;
    var yPos = options.yPos;
    var w = options.width;
    var h = options.height;
	var colourGradient = options.colourGradient;
	var fillPolygons = options.fillPolygons;
	var tooltips = options.tooltips;
	
	if (this.voronoi == undefined)
		this.voronoi = new greg.ross.visualisation.JSVoronoi(xPos, yPos, w, h, colourGradient, this.containerElement, fillPolygons, tooltips);
		
	this.voronoi.redraw(data);
}

/*
 * This class does most of the work.
 * *********************************
 */
greg.ross.visualisation.JSVoronoi = function(x, y, width, height, colourGradient, targetElement, fillRegions, tooltips)
{
	this.targetDiv;
	var id = allocateId();
	var canvas;
	var canvasContext = null;
	var minX;
	var minY;
	var maxX;
	var maxY;
	
	// The three intitial generators for the Voronoi diagram.
	var p1;
	var p2;
	var p3;
	
	// All generators.
	var generators = new Array();
	
	// Store vertices.
	var vertices = new Array();
	
	// Store polygons.
	var polygons = new Array();
	
	// Store edges.
	var edges = new Array();
	
	// Edges to be deleted.
	var delEdges = new Array();
	
	var qTree;
	var scaleX  = 1.9;  scaleY  = 1.9;
	var offsetX = -0.5;  offsetY = -0.5;
	
	function init()
    {
		createTargetDiv(); 
		
		if (!targetDiv) 
            return;
		
		createCanvas();
    }
	
	function render(data)
	{
		canvasContext.fillStyle = '#000'; 
		calculateMinMaxCoords(data);
		createSeedGraph();
		
		// Get the vertex in the centre of the circle that
		// circumbscribes the three points.
		createFirstVertex();
		
		// Create the vertices that represent the points on the
		// closed curve that surround the augmented geometric graph.
		createInfiniteVertices();
		
		// Create the winged-edge data structure.
		createWingedEdge();
		
		// Remove the temporary store of vertices.
		vertices = [];
		
		// Build the quaternary tree for fast indexing.
		
		var generators = new Array();
		var numPoints = data.getNumberOfRows();
		var i = numPoints - 1;
		var margin = 50;
		
		do 
	    {
	        var xCoord = data.getValue(i, 0);
			var yCoord = data.getValue(i, 1);
			xCoord = mapValueToZeroOneInterval(xCoord, minX, maxX);
			yCoord = mapValueToZeroOneInterval(yCoord, minY, maxY);
			
			var v = new greg.ross.visualisation.Vertex(xCoord, yCoord, 1);
			generators[i] = v;
	    }
	    while (i-- > 0)
		
		qTree = new greg.ross.visualisation.QuaternaryTree(generators, p1);
		
		try
		{
			populateDiagram();
		}
		catch(err)
		{
			//alert("Some edges are co-incident or parallel.");
		}
		
		if (fillRegions == true)
			fillPolygons();
		else	
			renderEdges();
			
		renderGenerators(data);
	}
	
	this.redraw = function(data)
	{
		render(data);
	}
	
	function fillPolygons()
	{
		canvasContext.lineWidth = 1;
		canvasContext.fillStyle = '#0000ff';
		canvasContext.strokeStyle='#888';
		var x1, x2, y1, y2;
		var n = polygons.length;
		
		for (var i = 0; i < n; i++)
		{
			var v = new Array();
			getEdgesAndVertices(polygons[i], null, v);
			canvasContext.beginPath();
			
			x1 = calcPosX(v[0].x);
			y1 = calcPosY(v[0].y);
			canvasContext.moveTo(x1, y1);
			var maxDistanceToGenerator = 0.0;
			var generator = (polygons[i]).generator;
			var nv = v.length;
			
			for (var j = 1; j < nv; j++)
			{
				if (v[j].w == 1)
				{
					x2 = calcPosX(v[j].x);
					y2 = calcPosY(v[j].y);
					
					canvasContext.lineTo(x2, y2);
					
					if (generator != null)
					{
						var dist = getDistance(v[j], generator);
						
						if (dist > maxDistanceToGenerator)
							maxDistanceToGenerator = dist;
					}
				}
			}
			
			if (generator != null)
			{
				var gx = calcPosX(generator.x);
				var gy = calcPosY(generator.y);
				var ieScaler = 1;
				
				if (isIE())
					ieScaler = 8;
					
				var radgrad = canvasContext.createRadialGradient(gx, gy, 0, gx, gy, calcPosX(Math.min(maxDistanceToGenerator, 0.2))*ieScaler);
					
				radgrad.addColorStop(0, colourGradient[0]);
				radgrad.addColorStop(0.9, colourGradient[1]);
				radgrad.addColorStop(1, colourGradient[2]);
				
				canvasContext.fillStyle = radgrad;
			}
			
			canvasContext.fill();
			canvasContext.closePath();
			canvasContext.stroke();
		}
	}
	
    function isIE()
    {
		return /msie/i.test(navigator.userAgent) && !/opera/i.test(navigator.userAgent);
    }
	
	/**
	* Render all edges in the Voronoi diagram
	*/
	
	function renderEdges()
	{
		var startV, endV;
		var x1, x2, y1, y2;
		canvasContext.beginPath();
		canvasContext.lineWidth = 1;
		canvasContext.strokeStyle='#888';
		canvasContext.lineJoin = "round";
		var n = edges.length;
		
		for (var i = 0; i < n; i++)
		{
			if (edges[i] != null)
			{
				if (!edges[i].deleted)
				{
					startV = edges[i].startVertex;
					endV = edges[i].endVertex;
					
					if ((startV.w == 1) && (endV.w == 1))
					{
						x1 = calcPosX(startV.x);
						y1 = calcPosY(startV.y);
						x2 = calcPosX(endV.x);
						y2 = calcPosY(endV.y);
						
						canvasContext.moveTo(x1, y1);
						canvasContext.lineTo(x2, y2);
					}
				}
			}
		}
		
		canvasContext.stroke();
	}
	
	/**
	* Render all generators (points) in the Voronoi diagram
	*/
	
	function renderGenerators(data)
	{
		// Render the generators.
		var numPoints = data.getNumberOfRows();
		var i = numPoints - 1;
		var margin = 50;
		canvasContext.fillStyle = '#ff2222';
		
		do 
	    {
	        var xCoord = data.getValue(i, 0);
			var yCoord = data.getValue(i, 1);
			
			xCoord = calcPosX(mapValueToZeroOneInterval(xCoord, minX, maxX));
			yCoord = calcPosY(mapValueToZeroOneInterval(yCoord, minY, maxY));
			
			canvasContext.beginPath();
			canvasContext.arc(xCoord, yCoord, 1, 0, self.Math.PI*2, true);
			canvasContext.fill();
	    }
	    while (i-- > 0)
	}
	
	/**
	* calculates the x position for this object within
	* the canvas and layout bounds of the model
	*
	*/
	
	function calcPosX(x)
	{
		return width / 2.0 + ((x + offsetX) * scaleX * (width / 2.0));
	}
	
	/**
	* calculates the y position for this object within
	* the canvas and layout bounds of the model
	*
	*/
	
	function calcPosY(y)
	{
		return height - (height/2.0 + ((y + offsetY) * scaleY * (height / 2.0)));
	}
	
	/**
	* Given a 2D point, determine the closest leaf node in the quaternary tree.
	*/
	function getClosestLeafNode(p, node)
	{
		var result = -1;
		var dist, closestDist = 1000000.0;
		var testNode = null;
		var closestNode = null;
		var scaledPoint = null;
		
		// Find the populated node which has a generator point closest to p
		if (node.nodeType != greg.ross.visualisation.QuadTreeNode.LEAF_NODE)
		{
			for (var j = 0; j < 4; j++)
			{
				// Determine which quadrant node we're examining
				
				switch(j)
				{
					case 0:
						testNode = node.topLeft;
						break;
					case 1:
						testNode = node.topRight;
						break;
					case 2:
						testNode = node.bottomLeft;
						break;
					case 3:
						testNode = node.bottomRight;
						break;
				}
				
				if (testNode.getGeneratorCount() > 0)
				{
					// Determine the closest node
					
					var tn = testNode.getGenerator(0);
					
					var scaledNodePoint = new greg.ross.visualisation.Vertex(calcPosX(tn.x), calcPosY(tn.y), 1);
					dist = getDistance(p, tn);
					
					if (dist < closestDist)
					{
						closestDist = dist;
						closestNode = testNode;
					}
				}
			}
			
			if (closestNode != null)
				result = getClosestLeafNode(p, closestNode);
			else
				result = getClosestLeafNode(p, testNode);
		}
		else
		{	
			result = node;
		}
		
		return result;
	}
	
	/**
	* Create the augmented geometric graph corresponding to the
	* Voronoi diagram of the intial three generator points.
	* These generator points contain the unit square containing the
	* rest of the Voronoi diagram.
	*/
	function createSeedGraph()
	{
		// Three points whose convex hull forms a triangle
		// which contains the unit square.
		
		p1 = new greg.ross.visualisation.Vertex(0.5, 3.0*Math.sqrt(2.0)/2.0 + 0.5, 1);
		
		p2 = new greg.ross.visualisation.Vertex(-3.0*Math.sqrt(6.0)/4.0 + 0.5, -3.0*Math.sqrt(2.0)/4.0 + 0.5, 1);
		
		p3 = new greg.ross.visualisation.Vertex(3.0*Math.sqrt(6.0)/4.0 + 0.5, -3.0*Math.sqrt(2.0)/4.0 + 0.5, 1);
		
		generators[0] = p1;
		generators[1] = p3;
		generators[2] = p2;
	}
	
	/**
	* From the centre point of the above circle, create
	* the first vertex of the Voronoi diagram.
	*/
	function createFirstVertex()
	{
		var centre = new greg.ross.visualisation.Vertex();
		
		if (getCenter(centre) == true)
		{
			var v = new greg.ross.visualisation.Vertex();
			v.x = centre.x;
			v.y = centre.y;
			v.w = 1;
			
			vertices[0] = v;
		}
	}
	
	/**
	* Given the three intial generators of the seed diagram,
	* determine the circle that circumbscribes them.
	* The centre of this circle is the coordinate of the first
	* non-infinite vertex.
	*/
	function getCenter(centre)
	{
		// Given three points, return the coordinates of the center of the circle passing
		// through them, return false if no such circle exists.
		var result = true;
		var x, y;
		
		var ma = 0;
		var mb = 0;
		
	    var nx1, nx2, nx3, ny1, ny2, ny3, ny4;
		
		nx1 = p1.x; 
		nx2 = p2.x; 
		nx3 = p3.x;
		ny1 = p1.y;
		ny2 = p2.y;
		ny3 = p3.y;
		
		if (nx1 != nx2)
			ma = (ny2 - ny1) / (nx2 - nx1);
		else
			result = false;
		
		if (nx2 != nx3)
			mb = (ny3 - ny3) / (nx3 - nx2);
		else
			result = false;
		
		if ((ma == 0) && (mb == 0))
			result = false;
		
		if (ma == mb)
			result = false;
		
		if (result == true)
		{
			x = (ma * mb * (ny1 - ny3) + mb * (nx1 + nx2) - ma * (nx2 + nx3)) / (2 * (mb - ma));
			
			if (ma != 0)
				y = -(x - (nx1 + nx2) / 2) / ma + (ny1 + ny2) / 2;
			else
				y = -(x - (nx2 + nx3) / 2) / mb + (ny2 + ny3) / 2;
			
			centre.x = x;
			centre.y = y;
		}
		
	    return result;
	}
	
	/**
	* Determine the infinite vertices on the closed curve containing the augmented
	* geometric graph. Do this by finding the line from the first non-infinite
	* vertex to the centre point of each line between the seed generators.
	*/
	function createInfiniteVertices()
	{
		addSeedVertex(p1, p3);
		addSeedVertex(p2, p3);
		addSeedVertex(p1, p2);
	}
	
	function addSeedVertex(c1, c2)
	{
		var midX, midY;
		var centreX = vertices[0].x;
		var centreY = vertices[0].y;
		var length;
		
		// Mid-point of line connecting c1 and c2.
		midX = (c1.x + c2.x) / 2.0;
		midY = (c1.y + c2.y) / 2.0;
		
		// Vector representing direction from centre vertex to this point.
		length = Math.sqrt(((midX - centreX) * (midX - centreX)) + ((midY - centreY) * (midY - centreY)));
		
		var v = new greg.ross.visualisation.Vertex();
		v.x = (midX - centreX) / length;
		v.y = (midY - centreY) / length;
		v.w = 0;
		
		vertices[vertices.length] = v;
	}
	
	/**
	* Given the vertices of the seed Voronoi diagram, derive the 
	* winged edge data structure.
	*/
	function createWingedEdge()
	{
		createEdges();
		createPolygons();
		createVertices();
		
		// Finish setting the polygon and vertex variables for the edges.
		edges[0].cwPredecessor = edges[2];
		edges[0].ccwPredecessor = edges[1];
		edges[0].cwSuccessor = edges[4];
		edges[0].ccwSuccessor = edges[3];
		edges[0].leftPolygon = polygons[1];
		edges[0].rightPolygon = polygons[0];
		
		edges[1].cwPredecessor = edges[0];
		edges[1].ccwPredecessor = edges[2];
		edges[1].cwSuccessor = edges[5];
		edges[1].ccwSuccessor = edges[4];
		edges[1].leftPolygon = polygons[2];
		edges[1].rightPolygon = polygons[1];
		
		edges[2].cwPredecessor = edges[1];
		edges[2].ccwPredecessor = edges[0];
		edges[2].cwSuccessor = edges[3];
		edges[2].ccwSuccessor = edges[5];
		edges[2].leftPolygon = polygons[0];
		edges[2].rightPolygon = polygons[2];
		
		edges[3].cwPredecessor = edges[5];
		edges[3].ccwPredecessor = edges[2];
		edges[3].cwSuccessor = edges[0];
		edges[3].ccwSuccessor = edges[4];
		edges[3].leftPolygon = polygons[0];
		edges[3].rightPolygon = polygons[3];
		
		edges[4].cwPredecessor = edges[3];
		edges[4].ccwPredecessor = edges[0];
		edges[4].cwSuccessor = edges[1];
		edges[4].ccwSuccessor = edges[5];
		edges[4].leftPolygon = polygons[1];
		edges[4].rightPolygon = polygons[3];
		
		edges[5].cwPredecessor = edges[4];
		edges[5].ccwPredecessor = edges[1];
		edges[5].cwSuccessor = edges[2];
		edges[5].ccwSuccessor = edges[3];
		edges[5].leftPolygon = polygons[2];
		edges[5].rightPolygon = polygons[3];
	}
	
	function populateDiagram()
	{
		var k = qTree.levelCount();
		var nodes;
		var node, parentNode;
		
		// Get each level of the tree
		
		for (var i = 0; i < k; i++)
		{
			nodes = qTree.getLevel(i);
			
			// If the node is not a leaf node, then only add the
			// generator to the diagram if it is not
			// the same generator that its parent references, because that's already
			// been added.
			var n = nodes.length;
			
			for (var j = 0; j < n; j++)
			{
				node = nodes[j];
				parentNode = node.parentNode;
				
				if (node.nodeType != greg.ross.visualisation.QuadTreeNode.LEAF_NODE)
				{
					if ((node.getGenerator(0) != null) &&
						node.getGenerator(0) != parentNode.getGenerator(0))
					{
						addGenerator(node.getGenerator(0), node);
					}
				}
				else
				{
					for (var m = 0; m < node.getGeneratorCount(); m++)
					{
						if (node.getGenerator(m) != parentNode.getGenerator(0))
							addGenerator(node.getGenerator(m), node);
					}
				}
			}
		}
	}
	
	/**
	* Given a coordinate, determine the closest generator to it
	*/
	function findClosestGenerator(p, node)
	{
		// Use quaternary initial guessing to find a close candidate
		// generator for starting the search for the actual closest one
		
		// If 'node' is not a leaf node then use its parent as the
		// initial guess
		
		var pi = null;
		
		if (node.nodeType != greg.ross.visualisation.QuadTreeNode.LEAF_NODE)
			pi = node.parentNode.getGenerator(0);
		else
			if (node.getGeneratorCount() > 1)
			{
				pi = node.getGenerator(0);
				
				if (pi.polygon == null)
					pi = node.parentNode.getGenerator(0);
			}
			else
				pi = node.parentNode.getGenerator(0);
		
		return getClosestPolygon(pi, p, false);
	}
	
	function getClosestPolygon(intitialGuess, p,  scalePoints)
	{
			var pk = null;
			var i;
			var pi, pj = null;
			
			var temp1;
			
			pi = intitialGuess;
			
			// Find the neighbouring generator pi that is closest
			// to p
			
			var e = new Array();
			var bFinished = false;
			
			while (!bFinished)
			{
				// Get edges around pi
				
				e = [];
				var pgi = pi.polygon;
				
				getEdgesAndVertices(pgi, e, null);
				
				// From each of these edges get the polygon that is
				// not pgi. Find the polygon that gives the smallest
				// distance to p
				
				var minDist = 999999.9;
				var dist;
				var otherP;
				var ed;
				var otherG;
				var ne = e.length;
				
				for (i = 0; i < ne; i++)
				{
					// Get the polygon
					
					ed = e[i];
					if (ed.leftPolygon == pgi)
						otherP = ed.rightPolygon;
					else
						otherP = ed.leftPolygon;
					
					otherG = otherP.generator;
					
					if (otherG != null)
					{
						pk = otherG;
						
						// Rescale the points so that they're in the coordinate
						// space of the view
						
						if (scalePoints)
						{
							temp1 = new greg.ross.visualisation.Vertex();
							temp1.x = (calcPosX(pk.x));
							temp1.y = (calcPosY(pk.y));
							temp1.w = pk.w;
						}
						else
							temp1 = pk;
						
						dist = getDistance(temp1, p);
						
						if (dist <= minDist)
						{
							minDist = dist;
							pj = otherG;
						}
					}
				}
				
				// Rescale the points so that they're in the coordinate
				// space of the view
				
				if (scalePoints)
				{
					temp1 = new greg.ross.visualisation.Vertex();
					temp1.x = (calcPosX(pi.x));
					temp1.y = (calcPosY(pi.y));
					temp1.w = pi.w;
				}
				else
					temp1 = pi;
				
				if (getDistance(temp1, p) <= minDist)
					bFinished = true;
				else
					pi = pj;
			}
				return  pi.polygon;
	}
	
	/**
	* Given a vertex j, determine the edges and polygons
	* that are incident to it. e is the list of edges incident to the vertex
	* and p is the list of polygons incident to the vertex
	*/
	function getEdgesAndPolygons(j, e, p)
	{
		var k, kStart;
		
		k = j.edgeAroundVertex;
		kStart = j.edgeAroundVertex;
		
		if (e != null)
			e[e.length] = k;
		
		var bFinished = false;
		
		while (bFinished == false)
		{
			if (j == k.startVertex)
			{
				if (p != null)
					p[p.length] = k.leftPolygon;
				
				k = k.ccwPredecessor;
			}
			else
			{
				if (p != null) 
					p[p.length] = k.rightPolygon;
				
				k = k.ccwSuccessor;
			}
			
			if (k == kStart)
			{
				bFinished = true;
			}
			else
				if (e != null)
					e[e.length] = k;
		}
	}
	
	/**
	* Given a polygon, determine the edges and vertices that
	* suround it. e is the list of edges around the polygon
	* and v is the list of vertices around the polygon
	*/
	function getEdgesAndVertices(p, e, v)
	{
		var k, kStart;
		
		k = p.edgeAroundPolygon;
		kStart = p.edgeAroundPolygon;
		
		if (e != null)
			e[e.length] = k;
		
		var bFinished = false;
				
		while (!bFinished)
		{
			if (p == k.leftPolygon)
			{
				if (v != null)
					v[v.length] = k.endVertex;
				
				k = k.cwSuccessor;
			}
			else
			{
				if (v != null)
					v[v.length] = k.startVertex;
				
				k = k.cwPredecessor;
			}
			
			if (k == kStart)
			{
				bFinished = true;
			}
			else
				if (e != null) 
					e[e.length] = k;
		}
	}
	
	/**
	* H(Pi, Pj, Pk, Pi) = 0 represents the circle passing through the points
	* Pi, Pj, Pk counterclockwise in this order. If H(Pi, Pj, Pk, P) < 0 then
	* the circle contains point P. If H(Pi, Pj, Pk, P) >= 0 where P is not in
	* {Pi, Pj, Pk}, then the circle is empty
	*/
	function H(Pi, Pj, Pk, Pl)
	{
		var J2ijk, J3ijk, J4ijk;
		var Xi, Xj, Xk, Xl, Yi, Yk, Yj, Yl;
		
		Xi = Pk.x;  Yk = Pi.y;
		Xj = Pj.x;  Yj = Pj.y;
		Xk = Pi.x;  Yi = Pk.y;
		Xl = Pl.x;  Yl = Pl.y;
		
		J2ijk = ((Yi - Yk) * (((Xj - Xk) * (Xj - Xk)) + ((Yj - Yk) * (Yj - Yk)))) - 
			((Yj - Yk) * (((Xi - Xk) * (Xi - Xk)) + ((Yi - Yk) * (Yi - Yk))));
			
		J3ijk = ((Xi - Xk) * (((Xj - Xk) * (Xj - Xk)) + ((Yj - Yk) * (Yj - Yk)))) - 
			((Xj - Xk) * (((Xi - Xk) * (Xi - Xk)) + ((Yi - Yk) * (Yi - Yk))));
			
		J4ijk = ((Xi - Xk) * (Yj - Yk)) - ((Xj - Xk) * (Yi - Yk));
		
		var result = ((J2ijk * (Xl - Xk)) - (J3ijk * (Yl - Yk)) + (J4ijk * (((Xl - Xk) * (Xl - Xk)) + ((Yl - Yk) * (Yl - Yk))))); 
		
		return result;
	}
	
	/**
	* Method to add a new generator to the Voronoi diagram
	*/
	function addGenerator(p, node)
	{
		// Get the polygon that surrounds the closest generator to p
		var closestG = findClosestGenerator(p, node);
		
		if (closestG != null)
		{
			var closestPolygon = closestG;
			
			// List of vertices surrounding the polygon
			var v = new Array();
			
			getEdgesAndVertices(closestPolygon, null, v);
			
			// Array to hold vertices that will be deleted
			var T = new Array();
			
			// For each vertex (Qijk) that is on the boundary of the polygon
			// find the one that gives the smallest value of H(Pi, Pj, Pk, P)
			// Add this vertex to T
			var minH = 999999.9;
			var hValue;
			var smallestV = null;
			var pg = new Array();
			var Pi, Pj, Pk;
			var nv = v.length;
			
			for (var i = 0; i < nv; i++)
			{
				pg = [];
				
				getEdgesAndPolygons(v[i], null, pg);
				
				Pi = pg[0].generator;
				Pj = pg[1].generator;
				Pk = pg[2].generator;
				
				// Make sure that we don't try to use the virtual generator
				// that's in the outermost, infinite region
				if ((Pi != null) && (Pj != null) && (Pk != null))
				{
					hValue = H(Pi, Pj, Pk, p);
					
					if (hValue < minH)
					{
						minH = hValue;
						smallestV = v[i];
					}
				}
			}
			
			T[T.length] = smallestV;
			
			// Expand the set of vertices T that will form the tree
			// that will be deleted to make way for the new Voronoi
			// polygon
			expandT(T, p, 0);
			deleteEdges();
			createNewVertices(T, closestPolygon, p);
			
			generators[generators.length] = p;
			
			T = [];
		}
	}
	
	/** for each vertex that is connected by an edge to the vertices in T,
	* add that vertex if H(Pi, Pj, Pk, P) < 0 and the vertex is not on the
	* outermost circuit of the geometric graph (v.getWi() = 0)
	*/
	function expandT(T, p, startNum)
	{
		var e;
		var startVertex, endVertex, vertexToAdd, vInT;
		
		vInT = T[startNum];
		
		var eg = new Array();
		
		getEdgesAndPolygons(vInT, eg, null);
		var neg = eg.length;
		
		for (var i = 0; i < neg; i++)
		{
			e = eg[i];
			
			if (e != null)
			{
				startVertex = e.startVertex;
				endVertex = e.endVertex;
				
				// Make sure that the edge is not wholly on the outermost circuit
				if ((startVertex.w == 1) || (endVertex.w == 1))
				{
					// Make sure that the start and end vertices are
					// not already in T
					if (!(contains(T, startVertex) && contains(T, endVertex)))
					{
						// The vertex should not be on the outermost circuit
						// of G
						if (vInT == startVertex)
							vertexToAdd = endVertex;
						else
							vertexToAdd = startVertex;
						
						if (vertexToAdd.w == 1)
						{
							// If H(Pi, Pj, Pk, P) < 0 then add
							// the vertex to T
							var pg = new Array();
							
							getEdgesAndPolygons(vertexToAdd, null, pg);
							var Pi, Pj, Pk;
							
							Pi = pg[0].generator;
							Pj = pg[1].generator;
							Pk = pg[2].generator;
							
							// Make sure that we don't try to use the virtual generator
							// that's in the outermost, infinite region
							if ((Pi != null) && (Pj != null) && (Pk != null))
							{
								if (H(Pi, Pj, Pk, p) < 0)
								{
									T[T.length] = vertexToAdd;
									expandT(T, p, (T.length - 1));
								}
							}
						}
						else
							vertexToAdd = null;
					}
					else
					{
						// Delete old edges that were completely contained by T
						delEdges[delEdges.length] = e;
					}
				}
			}
		}
	}
	
	/**
	* For every edge connecting a vertex in T with a vertex
	* not in T, create a new vertex on the edge and thus divide
	* the edge into two edges
	*/
	function createNewVertices(T, closestPolygon, p)
	{
		// The array list of edges that are to be split in two
		var modifiedEdges = new Array();
		var e;
		
		var i;
		
		var endV, startV;
		var ed = new Array();
		var nt = T.length;
		
		for (i = 0; i < nt; i++)
		{
			ed = [];
			getEdgesAndPolygons(T[i], ed, null);
			var ned = ed.length;
			
			for (var j = 0; j < ned; j++)
			{
				startV = ed[j].startVertex;
				endV = ed[j].endVertex;
				
				if (!(contains(T, startV) && contains(T, endV)))
					modifiedEdges[modifiedEdges.length] = ed[j];
			}
		}
		
		traceEdges(modifiedEdges, closestPolygon, p, T);
	}
	
	function traceEdges(modifiedEdges, closestPolygon, p, T)
	{
		var e, newEdge, first;
		
		// Create a placeholder for the new polygon
		var newPolygon = new greg.ross.visualisation.Polygon();
		newPolygon.generator = p;
		p.polygon = newPolygon;
		
		// Start with the edges that have the closest polygon to the new
		// generator, on their getLeftPolygon and/or getRighPolygon properties
		var firstEdge = null, secondEdge = null;
		var i, firstIndex = 0, secondIndex = 0;
		var n = modifiedEdges.length;
		
		for (i = 0; i < n; i++)
		{
			e = modifiedEdges[i];
			if (e.leftPolygon == closestPolygon)
			{
				if (firstEdge == null)
				{
					firstEdge = e;
					firstIndex = i;
				}
				else if (firstEdge != null)
				{
					secondEdge = e;
					secondIndex = i;
				}
				
				if ((firstEdge != null) && (secondEdge != null))
					break;
			}
			else if (e.rightPolygon == closestPolygon)
			{
				if (secondEdge == null)
				{
					secondEdge = e;
					secondIndex = i;
				}
				else if (secondEdge != null)
				{
					firstEdge = e;
					firstIndex = i;
				}
				
				if ((firstEdge != null) && (secondEdge != null))
					break;
			}
		}
		
		var newEdges = new Array(); // List of newly created edges
		
		// Create the first new edge
		newEdge = createFirstNewEdge(firstEdge, secondEdge, closestPolygon, p, newPolygon);
		
		if (newEdge == null)
			return;
		
		var aTemp = sortInitialEdges(newEdge, p, firstEdge, secondEdge);
		
		if (aTemp != null)
		{
			var temp = secondIndex;
			secondIndex = firstIndex;
			firstIndex = temp;
			
			firstEdge = aTemp[0];
			secondEdge = aTemp[1];
		}
		
		modifiedEdges[secondIndex] = null;
		
		first = firstEdge;
		
		newEdges[newEdges.length] = newEdge; // Store the new edge
		newPolygon.edgeAroundPolygon = newEdge; // Set the edge attribute for the new polygon
		
		// Create the other new edges
		
		var modEdges = new Array();
		modEdges[modEdges.length] = firstEdge;
		modEdges[modEdges.length] = secondEdge;
		
		var bFinished = false;
		var lastBounded = closestPolygon;
		var nextPolygon;
		
		while (!bFinished)
		{
			var prevStart = newEdges[newEdges.length - 1].endVertex;
			
			if (secondEdge.leftPolygon == lastBounded)
				nextPolygon = secondEdge.rightPolygon;
			else
				nextPolygon = secondEdge.leftPolygon;
				
			lastBounded = nextPolygon;
			var nm = modifiedEdges.length;
			
			for (i = 0; i < nm; i++)
			{
				// Search for the edge that shares a polygon other than the
				// the latest one bounded, with secondEdge
				
				e = modifiedEdges[i];
				
				if (e == null)
					continue;
				
				if ((e.leftPolygon == nextPolygon) || 
					(e.rightPolygon == nextPolygon))
				{
					secondIndex = i;
					firstEdge = secondEdge;
					secondEdge = e;
					break;
				}
			}
			
			modifiedEdges[secondIndex] = null;
			
			if (secondEdge == first)
			{
				// We've completed the loop bounding the subtree contained in T
				newEdge = new greg.ross.visualisation.Edge();
				newEdge.startVertex = prevStart;
				newEdge.endVertex = newEdges[0].startVertex;
				newEdge.leftPolygon = newPolygon;
				newEdge.rightPolygon = nextPolygon;
				newEdges[newEdges.length] = newEdge; // Store the new edge
				bFinished = true;
			}
			else
			{
				newEdge = createNewEdge(secondEdge,  nextPolygon, p, newPolygon);
				newEdge.startVertex = prevStart;
				newEdges[newEdges.length] = newEdge; // Store the new edge
				modEdges[modEdges.length] = secondEdge;
			}
		}
		
		// Add the new edges to the set of existing edges
		for (i = 0; i < newEdges.length; i++)
		{
			e = newEdges[i];
			edges[edges.length] = e;
		}
		
		modifyEdges(newEdges, modEdges, T, newPolygon);
		
		// Add the new polygon to the last entry of the polygons array
		
		polygons[polygons.length] = newPolygon;
	}
	
	/**
	* Given the set of edges that are to be split in half and
	* the newly created edges, split the modified edges and
	* delete the sub-tree contained in T
	*/
	function modifyEdges(newEdges, modEdges, T, newPolygon)
	{
		var eMod, eNew, e;
		
		// Modify existing edges
		var i;
		for (i = 0; i < modEdges.length; i ++)
		{
			eNew = newEdges[i];
			eMod = modEdges[i];
			
			if (contains(T, eMod.startVertex))
			{
				eMod.startVertex = eNew.startVertex;
				
				eMod.ccwPredecessor = eNew;
				
				if (i > 0)
					eMod.cwPredecessor = newEdges[i - 1];
				else
					eMod.cwPredecessor = newEdges[newEdges.length - 1];
			}
			else if (contains(T, eMod.endVertex))
			{
				eMod.endVertex = eNew.startVertex;
				
				eMod.ccwSuccessor = eNew;
				
				if (i > 0)
					eMod.cwSuccessor = newEdges[i - 1];
				else
					eMod.cwSuccessor = newEdges[newEdges.length - 1];
			}
		}
		
		// Update the predecessor and successor properties of the newly
		// created edges
		var nen = newEdges.length;
		
		for (i = 0; i < nen; i++)
		{
			e = newEdges[i];
			e.cwPredecessor = modEdges[i];
			
			if (i > 0)
				e.ccwPredecessor = newEdges[i - 1];
			else
				e.ccwPredecessor = newEdges[newEdges.length - 1];
			
			if (i == (newEdges.length - 1))
			{
				e.ccwSuccessor = modEdges[0];
				e.cwSuccessor = newEdges[0];
			}
			else
			{
				e.ccwSuccessor = modEdges[i + 1];
				e.cwSuccessor = newEdges[i + 1];
			}
		}
	}
	
	/**
	* Project a line between the generator of poli and p and then find
	* the coordinates of where the perpendicular bisector of this line
	* intersects firstEdge. These will be the coordinates of
	* the new vertice. Draw an edge between these
	*/
	function createNewEdge(secondEdge, poli,  p, 
	 newPolygon)
	{
		// Get the generator for polygon pi
		var pi = poli.generator;
		
		// Get the line between pi and p
		var x1, x2, y1, y2;
		x1 = p.x; x2 = pi.x;
		y1 = p.y;  y2 = pi.y;
		
		// Create a line representing the perpendicular bisector of this line
		var mid = new greg.ross.visualisation.Vertex((x1 + x2) / 2.0, (y1 + y2) / 2.0, 0); // Mid point of line
		var bi = new greg.ross.visualisation.Vertex((y1 - y2), -(x1 - x2), 0); // Direction of perpendicular line
		
		var px1, px2, py1, py2; // Segment of the bisector
		px1 = mid.x;
		py1 = mid.y;
		px2 = mid.x + bi.x;
		py2 = mid.y + bi.y;
		
		// Find the point where this line intersects firstEdge and secondEdge
		// then create the edge with these coordinates as the start and end
		// vertices
		var newEdge = new greg.ross.visualisation.Edge();
		var v;
		
		v = getIntersection(secondEdge, px1, py1, px2, py2);
		v.w = 1;
		
		newEdge.endVertex = v;
		newEdge.endVertex.edgeAroundVertex = newEdge;
		
		// Set the polygon attributes
		newEdge.rightPolygon = poli;
		newEdge.leftPolygon = newPolygon;
		
		return newEdge;
	}
	
	/**
	* When creating the first new Voronoi edge, we must determine
	* upon which of the two existing edges the new edge's startVertex will be.
	*/
	
	function sortInitialEdges(newEdge, p, firstEdge,
	 secondEdge)
	{
		// Make sure that p is to the left of the new edge
		var v1, v2;
		v1 = newEdge.startVertex;
		v2 = newEdge.endVertex;
		
		var result = null;
		
		var iTest = (p.y - v1.y) * (v2.x - v1.x) - 
			(p.x - v1.x) * (v2.y - v1.y);
		
		if (iTest > 0)
		{
			result = new Array();
			newEdge.startVertex = v2;
			newEdge.endVertex = v1;
			result[result.length] = secondEdge;
			result[result.length] = firstEdge;
		}
		
		return result;
	}
	
	/**
	* Project a line between the generator of poli and p and then find
	* the coordinates of where the perpendicular bisector of this line
	* intersects firstEdge and secondEdge. These will be the coordinates of
	* the new vertices. Draw an edge between these
	*/
	function createFirstNewEdge(firstEdge, secondEdge, poli,  p, 
	 newPolygon)
	{
		// Get the generator for polygon pi
		var pi = poli.generator;
		
		// Get the line between pi and p
		var x1, x2, y1, y2;
		x1 = p.x; x2 = pi.x;
		y1 = p.y;  y2 = pi.y;
		
		// Create a line representing the perpendicular bisector of this line
		var mid = new greg.ross.visualisation.Vertex((x1 + x2) / 2.0, (y1 + y2) / 2.0, 0); // Mid point of line
		var bi = new greg.ross.visualisation.Vertex((y1 - y2), -(x1 - x2), 0); // Direction of perpendicular line
		
		var px1, px2, py1, py2; // Segment of the bisector
		px1 = mid.x;
		py1 = mid.y;
		px2 = mid.x + bi.x;
		py2 = mid.y + bi.y;
		
		// Find the point where this line intersects firstEdge and secondEdge
		// then create the edge with these coordinates as the start and end
		// vertices
		var newEdge = new greg.ross.visualisation.Edge();
		
		var vTemp = getIntersection(firstEdge, px1, py1, px2, py2);
		
		if (vTemp == null)
			return null;
		
		var v = new greg.ross.visualisation.Vertex(vTemp.x, vTemp.y, 1);
		
		newEdge.startVertex = v;
		newEdge.startVertex.edgeAroundVertex = newEdge;
		
		var vTemp = getIntersection(secondEdge, px1, py1, px2, py2);
		
		if (vTemp == null)
			return null;
		
		var v2 = new greg.ross.visualisation.Vertex(vTemp.x, vTemp.y, 1);
		
		newEdge.endVertex = v2;
		newEdge.endVertex.edgeAroundVertex = newEdge;
		
		// Set the polygon attributes
		newEdge.rightPolygon = poli;
		newEdge.leftPolygon = newPolygon;
		
		return newEdge;
	}
	
	/**
	* Given an edge and a line, return a coordinate representing the point
	* at which the line intersects the edge
	*/
	function getIntersection(e, xa, ya, xb, yb)
	{
		var x1 = xa; var y1 = ya;
		var x2 = xb; var y2 = yb;
		var x3, x4, y3, y4;
		
		if (e.startVertex.w == 0)
		{
			x3 = e.endVertex.x + e.startVertex.x;
			x4 = e.endVertex.x;
			y3 = e.endVertex.y + e.startVertex.y;
			y4 = e.endVertex.y;
		}
		else if (e.endVertex.w == 0)
		{
			x3 = e.startVertex.x;
			x4 = e.startVertex.x + e.endVertex.x;
			y3 = e.startVertex.y;
			y4 = e.startVertex.y + e.endVertex.y;
		}
		else
		{
			x3 = e.startVertex.x;
			x4 = e.endVertex.x;
			y3 = e.startVertex.y;
			y4 = e.endVertex.y;
		}
		
		var nma, nmb;  	// The two numerators
		var dm;      	// The common denominator
		var ua, ub;    	// The coefficients.
		
		nma = (((x4 - x3) * (y1 - y3)) - ((y4 - y3) * (x1 - x3)));
		nmb = (((x2 - x1) * (y1 - y3)) - ((y2 - y1) * (x1 - x3)));
		dm = ((y4 - y3) * (x2 - x1)) - ((x4 - x3) * (y2 - y1));
		
		if ((dm == 0) && (nma == 0)){}
			//alert("Some Voronoi edges are co-incident.");
		else if (dm == 0){}
			//alert("Some Voronoi edges are parallel.");
		else
		{
			ua = nma / dm;
			ub = nmb / dm;
			
			var x, y;
			
			x = x1 + (ua * (x2 - x1));
			y = y1 + (ua * (y2 - y1));
			
			return new greg.ross.visualisation.Vertex(x, y, 0);
		}
		
		return null;
	}
	
	function contains(T, value)
	{
		var n = T.length;
		
		for (var i = 0; i < n; i++)
		{
			if (T[i] == value)
				return true;
		}
	}
	
	function deleteEdges()
	{
		var index;
		var  e;
		var n = delEdges.length;
		
		for (var i = 0; i < n; i++)
		{
			e = delEdges[i];
			
			updatePolygonForEdgeDeletion(e, e.leftPolygon, delEdges);
			updatePolygonForEdgeDeletion(e, e.rightPolygon, delEdges);
			
			e.deleted = true;
		}
		
		delEdges = [];
	}
	
	/** 
	* If the edge to be deleted is either of its
	* neighbouring polygons' edgeAroundPolygon attribute
	* then get the edges around that polygon and set
	* reset this attribute for another edge
	*/
	
	function updatePolygonForEdgeDeletion(delEdge, p, delEdges)
	{
		var ed = new Array();
		var e;
		var i;
		
		if (p.edgeAroundPolygon == delEdge)
		{
			getEdgesAndVertices(p, ed, null);
			var ned = ed.length;
			
			for (i = 0; i < ned; i++)
			{
				e = ed[i];
				if ((e != delEdge) && (!contains(delEdges, e)))
				{
					p.edgeAroundPolygon = e;
					break;
				}
			}
		}
	}
	
	function createVertices()
	{
		var v;
		
		v = vertices[0];
		v.edgeAroundVertex = edges[0];
		
		v = vertices[1];
		v.edgeAroundVertex = edges[3];
		
		v = vertices[2];
		v.edgeAroundVertex = edges[4];
		
		v = vertices[3];
		v.edgeAroundVertex = edges[5];
	}
	
	function createPolygons()
	{
		var p;
		
		// Initially there are four polygons, including
		// the theoretical outermost infinite region.
		for (var i = 0; i < 4; i++)
		{
			p = new greg.ross.visualisation.Polygon();
			p.edgeAroundPolygon = edges[i];
			polygons[polygons.length] = p;
		}
		
		polygons[0].generator = p1;
		polygons[1].generator = p3;
		polygons[2].generator = p2;
		
		p1.polygon = polygons[0];
		p3.polygon = polygons[1];
		p2.polygon = polygons[2];
	}
	
	function createEdges()
	{
		var e;
		
		// Create edge 0.
		e = new greg.ross.visualisation.Edge();
		e.startVertex = vertices[0];
		e.endVertex = vertices[1];
		edges[edges.length] = e;
		
		// Create edge 1
		
		e = new greg.ross.visualisation.Edge();
		e.startVertex = vertices[0];
		e.endVertex = vertices[2];
		edges[edges.length] = e;
		
		// Create edge 2
		
		e = new greg.ross.visualisation.Edge();
		e.startVertex = vertices[0];
		e.endVertex = vertices[3];
		edges[edges.length] = e;
		
		// Create edge 3
		
		e = new greg.ross.visualisation.Edge();
		e.startVertex = vertices[3];
		e.endVertex = vertices[1];
		edges[edges.length] = e;
		
		// Create edge 4
		
		e = new greg.ross.visualisation.Edge();
		e.startVertex = vertices[1];
		e.endVertex = vertices[2];
		edges[edges.length] = e;
		
		// Create edge 5
		
		e = new greg.ross.visualisation.Edge();
		e.startVertex = vertices[2];
		e.endVertex = vertices[3];
		edges[edges.length] = e;
	}
	
	function calculateMinMaxCoords(data)
	{
		var numPoints = data.getNumberOfRows();
		var i = numPoints - 1;
		
		do 
	    {
	        var xCoord = data.getValue(i, 0);
			var yCoord = data.getValue(i, 1);
			
			if (minX == undefined) minX = xCoord;
			if (minY == undefined) minY = yCoord;
			if (maxX == undefined) maxX = xCoord;
			if (maxY == undefined) maxY = yCoord;
			
			if (xCoord < minX) minX = xCoord;
			if (xCoord > maxX) maxX = xCoord;
			if (yCoord < minY) minY = yCoord;
			if (yCoord > maxY) maxY = yCoord;
	    }
	    while (i-- > 0)
	}
    
    function createTargetDiv()
	{
		this.targetDiv = document.createElement("div");
		this.targetDiv.id = id;
		this.targetDiv.className = "voronoi";
		this.targetDiv.style.background = '#ffffff'
		this.targetDiv.style.position = 'absolute';
		
		if (!targetElement) 
			document.body.appendChild(this.targetDiv);
		else 
		{
			this.targetDiv.style.position = 'relative';
			targetElement.appendChild(this.targetDiv);
		}
		
		this.targetDiv.style.left = x + "px";
		this.targetDiv.style.top = y + "px";
	}
	
	function createCanvas()
	{
		canvas = document.createElement("canvas");
		
		if (isIE())
		{
			G_vmlCanvasManager.initElement(canvas);
			canvas.style.width = width;
			canvas.style.height = height;
		}
			
		canvas.className = "voronoiCanvas";
        canvas.setAttribute("width", width);
        canvas.setAttribute("height", height);
		canvas.style.left = '0px';
		canvas.style.top =  '0px';
		
		targetDiv.appendChild(canvas);
		
		canvasContext = canvas.getContext("2d");
		canvasContext.clearRect(0, 0, canvas.width, canvas.height);
		
		if (fillRegions == true)
			canvasContext.fillStyle = colourGradient[1];
		else
			canvasContext.fillStyle = '#fff';
		
		canvasContext.fillRect(0, 0, canvas.width, canvas.height);
		
		canvasContext.beginPath();
		canvasContext.rect(0, 0, canvas.width, canvas.height);
		canvasContext.strokeStyle='#888';
		canvasContext.stroke();
		
		canvas.onmousemove = displayTooltip;
		canvas.onmouseout = hideTooltip;
	}
	
	function hideTooltip()
	{
		tooltip.hide();
	}
	
	function displayTooltip(e)
	{
		if (isIE())
			var e = window.event;

		var xi = (e.clientX - (width / 2.0)) / (scaleX * (width / 2.0)) - offsetX;
		var yi = ((height - e.clientY) - (height / 2.0)) / (scaleY * (height / 2.0)) - offsetY;
		
		var position = new greg.ross.visualisation.Vertex(e.clientX-x, e.clientY-y, 1);
		var rootNode = qTree.getRootNode();
		var closestLeafNode = getClosestLeafNode(position, rootNode);
		
		var closestPolygon = getClosestPolygon(closestLeafNode.getGenerator(0), position,  true)
		
		tooltip.show(tooltips[closestPolygon.generator.index], 200);
	}
	
	function allocateId()
	{
		var count = 0;
		var name = "voronoi";
		
		do
		{
			count++;
		}
		while(document.getElementById(name+count))
		
		return name+count;
	}
	
	function mapValueToZeroOneInterval(value, minValue, maxValue)
	{
		if (minValue == maxValue) return 0;
		
		var factor = (value - minValue) / (maxValue - minValue);
		return factor;
	}
	
	init();
}

/*
 * Vertex: represents a 2-d point on the plane.
 * When it is the defining point of a Voronoi region it is often
 * referred to as a "generator" in the nomenclature.
 * 
 * w: Contains 1 if vertex i is an ordinary point, otherwise
 * contains 0 if point i is a point at infinity
 *	
 * x and y: If Vertex i is an ordinary point,
 * that is w = 1 then store the x and y co-ordinates
 * respectively. On the other hand if w = 0 then store the
 * x and y components of the unit vector designating the direction
 * in which the associated (infinite) Voronoi edge runs
 * ******************************************************************
 */
greg.ross.visualisation.Vertex = function(x, y, w)
{
	this.x = x;
	this.y = y;
	this.w = w;
	this.index = 0;
	this.polygon = null;
	
	// edgeAroundVertex: An edge incident to the
	// vertex in the "winged edge data structure".
	this.edgeAroundVertex;
}

/*
 * Edge: As defined in the winged edge data structure.
 * ***************************************************
 */
greg.ross.visualisation.Edge = function()
{
	// rightPolygon: The polygon that is to
	// the right of the edge.
	this.rightPolygon;
	
	// leftPolygon: The polygon that is to
	// the left of the edge.
	this.leftPolygon;
	
	// startVertex: The start vertex of the
	// edge
	this.startVertex;
	
	// endVertex: The end vertex of the
	// edge
	this.endVertex;
	
	// cwPredecessor: The edge next to this edge
	// clockwise around the start vertex.
	this.cwPredecessor;
	
	// ccwPredecessor: The edge next to this edge
	// counterclockwise around the start vertex.
	this.ccwPredecessor;
	
	// cwSuccessor: The next to this edge
	// clockwise around the end vertex.
	this.cwSuccessor;
	
	// ccwSuccessor: The edge next to this edge
	// counterclockwise around the end vertex.
	this.ccwSuccessor;
	
	// Determine whether this edge should be deleted.
	this.deleted = false;
}

/*
 * Polygon: As defined in the winged edge data structure.
 * ***************************************************
 */
greg.ross.visualisation.Polygon = function()
{
	// edgeAroundPolygon: An edge on the boundary of
	// the polygon.
	this.edgeAroundPolygon;
	
	// The vertex that is associated with the polygon.
	this.generator = null;
}


greg.ross.visualisation.QuaternaryTree = function(generators, p1)
{
	// The number of levels in the tree
	var k = 0;
	
	// The number of generators
	var numPoints = 0;
	
	// Number of buckets (leaf nodes) in the tree
	var numBuckets = 0;
	
	// Array of buckets
	var buckets = new Array();
	
	// Array of tree levels. Each element is an array list of
	// tree nodes on the corresponding level
	var treeLevels = new Array();
	
	var rootNode;
	
	// The number of rows and columns in the grid representing
	// the entire set of leaf node
	var numRows;
	
	// The coordinate that will be assigned to the rootnode
	var p1;
	
	createTree();
	
	/**
	* Build the quaternary tree
	*/
	function createTree()
	{
		numPoints = generators.length;
		defineK();
		addPointsToBuckets();
		
		// Create the root node
		rootNode = new greg.ross.visualisation.QuadTreeNode();
		rootNode.nodeType = greg.ross.visualisation.QuadTreeNode.ROOT_NODE;
		rootNode.addGenerator(p1);
		
		// Build the quaternary tree
		buildTree(k, rootNode, 0, numRows - 1, 0, numRows - 1);
	}
	
	/**
	* Given the number of genereator points, determine the
	* number of levels and leaf nodes in the tree
	*/
	function defineK()
	{
		k = Math.round(Math.log(numPoints) / Math.log(4));
		numBuckets = Math.pow(4, k);
		numRows = Math.sqrt(numBuckets);
		
		// Init the array of tree levels
		treeLevels = new Array();
		
		for (var i = 0; i < k; i++)
			treeLevels[i] = new Array();
	}
	
	/**
	* Allocate each generator to a bucket in the tree
	*/
	function addPointsToBuckets()
	{
		buckets = new Array();
		var r, c, index;
		var k2 = Math.pow(2, k);
		var n = generators.length;
		
		for (var i = 0; i < n; i++)
		{
			// Determine the row and column by multiplying the
			// coordinate by 2^k and then truncating off the fractional
			// parts
			var g = generators[i];
			
			// Store a reference to the data point in the Vertex object.
			// This will facilitate view coordination.
			g.index = i;
			
			var x = 0, y = 0;
			if (g.x == 1) 
				x = 0.999999;
			else
				x = g.x;
			
			if (g.y == 1)
				y = 0.999999;
			else
				y = g.y;
			
			r = parseInt(k2 * y);
			c = parseInt(k2 * x);
			
			// The generators are added to a 1D array of buckets. Each
			// index position is derived from the r and c values
			var index = c + (r * numRows);
			
			if (buckets[index] == undefined)
				buckets[index] = new Array();
			
			(buckets[index])[buckets[index].length] = generators[i];
		}
	}
	
	/**
	* Build the structure of the tree in a depth-first manner
	*/
	function buildTree(kIn, n, rStart, rEnd, cStart, cEnd)
	{
		if (kIn > 0)
		{
			var newNode;
			var divider1 = parseInt(rStart + ((rEnd - rStart) / 2.0));
			var divider2 = parseInt(cStart + ((cEnd - cStart) / 2.0));
			
			for (var i = 0; i < 4; i++)
			{
				newNode = new greg.ross.visualisation.QuadTreeNode();
				
				switch(i)
				{
					case 0:
						n.topLeft = newNode;
						newNode.parentNode = n;
						buildTree(kIn - 1, newNode, rStart, divider1, cStart, divider2);
						break;
					case 1:
						n.topRight = newNode;
						newNode.parentNode = n;
						buildTree(kIn - 1, newNode, rStart, divider1, divider2 + 1, cEnd);
						break;
					case 2:
						n.bottomLeft = newNode;
						newNode.parentNode = n;
						buildTree(kIn - 1, newNode, divider1 + 1, rEnd, cStart, divider2);
						break;
					case 3:
						n.bottomRight = newNode;
						newNode.parentNode = n;
						buildTree(kIn - 1, newNode, divider1 + 1, rEnd, divider2 + 1, cEnd);
						break;
				}
				
				// Store the nodes for this particular level
				var idx = Math.abs(kIn - k);
				(treeLevels[idx])[treeLevels[idx].length] = newNode;
			}
		}
		else
		{
			// Add the corresponding buckets to each of the leaf nodes
			n.nodeType = greg.ross.visualisation.QuadTreeNode.LEAF_NODE;
			
			var index = parseInt(cStart + (rStart * numRows));
			
			n.setGenerators(buckets[index]);
			
			// Populate the parent nodes up to but not including the root node
			if (buckets[index] != null)
				populateParentNodes(n, (buckets[index])[0]);
		}
	}
	
	/**
	* After creating and populating the leaf nodes, we must select a
	* generator from each non-empty leaf node bucket and allocate it
	* to each consecutive parent defined on the path from the leaf node
	* to the root node
	*/
	
	function populateParentNodes(child, g)
	{
		var parent = child.parentNode;
		
		if ((parent.nodeType != greg.ross.visualisation.QuadTreeNode.ROOT_NODE) && (parent.getGeneratorCount() == 0))
		{
			parent.addGenerator(g);
			populateParentNodes(parent, g);
		}
	}
	
	/**
	* Acccessor method to get the number of levels in the tree
	*/
	this.levelCount = function()
	{
		return k;
	}
	
	/**
	* Accessor method to get nodes in a particular level of the tree
	*/
	this.getLevel = function(i)
	{
		if ((i < 0) || (i > (k - 1)))
			return null;
		else
			return treeLevels[i];
	}
	
	this.getRootNode = function()
	{
		return rootNode;
	}
}


greg.ross.visualisation.QuadTreeNode = function()
{
	// Pointer to the child nodes defined by quadrant. These are null for
	// leaf nodes
	this.topLeft = null;
	this.topRight = null;
	this.bottomLeft = null;
	this.bottomRight = null;
	
	// Pointer to the parent node. This is null for
	// the root node
	this.parentNode = null;
	
	// Generators associated with this node
	this.generators = null;
	
	this.nodeType;
	
	// Accessor methods for the generator(s) associated with this node
	this.getGenerator = function(index)
	{
		if (this.generators == null)
			return null;
		else
			return this.generators[index];
	}
	
	this.addGenerator = function(generator)
	{
		if (this.generators == null)
			this.generators = new Array();
		
		this.generators[this.generators.length] = generator;
	}
	
	this.setGenerators = function(g)
	{
		this.generators = g;
	}
	
	this.getGeneratorCount = function()
	{
		if (this.generators == null)
			return 0;
		else
			return this.generators.length;
	}
}

/**
* Given two coordinates, return the Euclidean distance
* between them
*/
function getDistance(p1, p2)
{
	return Math.sqrt(((p1.x - p2.x) * (p1.x - 
		p2.x)) + ((p1.y - p2.y) * (p1.y - p2.y)));
}

/*
 * This function displays tooltips and was adapted from original code by Michael Leigeber.
 * See http://www.leigeber.com/
 */
var tooltip = function()
{
	var id = 'tt';
	var top = 3;
	var left = 3;
	var maxw = 300;
	var speed = 10;
	var timer = 20;
	var endalpha = 95;
	var alpha = 0;
	var tt,t,c,b,h;
	var ie = document.all ? true : false;
	
	return{
		show:function(v,w)
		{
			if (tt == null)
			{
				tt = document.createElement('div');
				tt.setAttribute('id',id);
				
				tt.style.position = 'absolute';
				tt.style.display =  'block';
				
				t = document.createElement('div');
				t.setAttribute('id',id + 'top');
				
				t.style.display = 'block';
				t.style.height =  '5px';
				t.style.marginleft =  '5px';
				t.style.overflow =  'hidden';
				
				c = document.createElement('div');
				c.setAttribute('id',id + 'cont');
				
				b = document.createElement('div');
				b.setAttribute('id',id + 'bot');
				
				tt.appendChild(t);
				tt.appendChild(c);
				tt.appendChild(b);
				document.body.appendChild(tt);
				
				if (!ie)
				{
					tt.style.opacity = 0;
					tt.style.filter = 'alpha(opacity=0)';
				}
				else
					tt.style.opacity = 1;
				
				document.onmousemove = this.pos;
			}
			
			tt.style.display = 'block';
			c.innerHTML = '<span style="font-weight:bold; font-family: arial;">' + v + '</span>';
			tt.style.width = w ? w + 'px' : 'auto';
			
			if (!w && ie)
			{
				t.style.display = 'none';
				b.style.display = 'none';
				tt.style.width = tt.offsetWidth;
				t.style.display = 'block';
				b.style.display = 'block';
			}
			
			if (tt.offsetWidth > maxw)
			{
				tt.style.width = maxw + 'px';
			}
			
			h = parseInt(tt.offsetHeight) + top;
			
			if (!ie)
			{
				clearInterval(tt.timer);
				tt.timer = setInterval(function(){tooltip.fade(1)},timer);
			}
		},
		pos:function(e)
		{
			var u = ie ? event.clientY + document.documentElement.scrollTop : e.pageY;
			var l = ie ? event.clientX + document.documentElement.scrollLeft : e.pageX;
			tt.style.top = (u - h) + 'px';
			tt.style.left = (l + left) + 'px';
		},
		fade:function(d)
		{
			var a = alpha;
			
			if ((a != endalpha && d == 1) || (a != 0 && d == -1))
			{
				var i = speed;
				
				if (endalpha - a < speed && d == 1)
				{
					i = endalpha - a;
				}
				else if	(alpha < speed && d == -1)
				{
					i = a;
				}
			
			alpha = a + (i * d);
			tt.style.opacity = alpha * .01;
			tt.style.filter = 'alpha(opacity=' + alpha + ')';
			}
			else
			{
				clearInterval(tt.timer);
				
				if (d == -1)
				{
					tt.style.display = 'none';
				}
			}
		},
		hide:function()
		{
			if (tt == null)
				return;
		
			if (!ie)
			{
				clearInterval(tt.timer);
				tt.timer = setInterval(function(){tooltip.fade(-1)},timer);
			}
			else
			{
				tt.style.display = 'none';
			}
		}
	};
}();


greg.ross.visualisation.QuadTreeNode.ROOT_NODE = 0;
greg.ross.visualisation.QuadTreeNode.INTERMEDDIATE_NODE = 1;
greg.ross.visualisation.QuadTreeNode.LEAF_NODE = 2;