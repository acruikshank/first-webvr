/*
Step: [0-1]
StepSink: Step ->
StepSource: StepSink ->
  step(n): stepSource // produce n linear steps between 0 and 1 inclusive
VertexSink: Vertex, Step, Step* ->
VertexGenerator: VertexSink -> StepSink
  vertices(vertices) // take the given vertices to create a VertexGenerator
  parametric(f): (Step -> Vertex) -> VertexGenerator // convert a vertex producer into a vertex generator
VertexMap: VertexSink -> VertexSink
  zTranslate(start,end): translate vertices along z axis
  zRotate(start*,end*): rotate vertices around z axis
  tag(edges): mark vertices with tags
Facer: FaceSink -> VertexSink
  face(options): Produce a facer with the given options for tesselating and joining
Renderer: FaceSink
  ThreeJSRenderer
  CSGRenerer
  STLRenderer
  NativeRenderer
*/

(function(context) {
  var debug = false;
  /*
  // cube
  var renderer = ThreeJSRenderer();

  // cube
  step(2)(
    vertices([[1,1,0],[-1,1,0],[-1,-1,0],[1,-1,0]])(
    transform(zTranslate(1), label) (
    face( edgeLoop, tesselate('bottom'), tesselate('top',true) )( renderer ))))

  var geometry = renderer.geometry;

  // OR
  extrudeClosed = face( edgeLoop, tesselate('bottom'), tesselate('top',true) )

  // sphere
  var R = 40, Q = 20, 2PI = 2*Math.Pi, hPI = Math.PI/2;
  function simiCircle(s) {return [R*Math.cos(tween(s,-hPI,hPI)),0,R*Math.sin(tween(s,-hPI,hPI))]}

  step(Q)(
    step(Q)(
      parametric(semiCircle))(
        zRotate()(
          face(connect:'tb', singular:'lr')( renderer ))))



  */
  var vertexId = 1;

  function Vertex( point, transformStep, ribStep, id) {
    this[0] = point[0];
    this[1] = point[1];
    this[2] = point[2];
    this[3] = transformStep;
    this[4] = ribStep;
    this[5] = id == null ? ++vertexId : id;
  }
  Vertex.prototype = new Float64Array(6);
  Object.defineProperty(Vertex.prototype, 'x', {get:function() {return this[0]}, set:function(x) {this[0]=x}})
  Object.defineProperty(Vertex.prototype, 'y', {get:function() {return this[1]}, set:function(y) {this[1]=y}})
  Object.defineProperty(Vertex.prototype, 'z', {get:function() {return this[2]}, set:function(z) {this[2]=z}})
  Object.defineProperty(Vertex.prototype, 'transformStep', {
    get:function() {return this[3]},
    set:function(transformStep) {this[3]=transformStep}})
  Object.defineProperty(Vertex.prototype, 'ribStep', {
    get:function() {return this[4]},
    set:function(ribStep) {this[4]=ribStep}})
  Object.defineProperty(Vertex.prototype, 'id', {
    get:function() {return this[5]},
    set:function(id) {this[5]=id}})
  Vertex.prototype.toString = function() { return "Vertex("+this.id+") ["+this.x+","+this.y+","+this.z+"] "+
                                                  this.transformStep+" "+this.ribStep };


  function resolveCurve( points, s, start, size ) {
    start = start || 0;
    size = size || points.length;
    if ( isNaN(size) || size < 2 ) return points[start];

    var p1 = resolveCurve( points, s, start, size-1 ), p2 = resolveCurve( points, s, start+1, size-1 );
    return vadd( p1, vscale( vsub(p2, p1), s) );
  }

  function curveTangent( points, s, start, size ) {
    start = start || 0;
    size = size || points.length;
    if ( size < 2 ) return [0,0,0];

    var p1 = resolveCurve( points, s, start, size-1 ), p2 = resolveCurve( points, s, start+1, size-1 );
    return vnorm( vsub( p2, p1 ) );
  }

  function curveCurl( points, s, start, size ) {
    start = start || 0;
    size = size || points.length;
    if ( size < 3 ) return [0,0,0];

    var p1 = resolveCurve( points, s, start, size-2 ),
        p2 = resolveCurve( points, s, start+1, size-2 ),
        p3 = resolveCurve( points, s, start+2, size-2 );
    return vcross( vsub( p3, p2 ), vsub( p1, p2 ) );
  }

  function splitCurve( points, s, start, size ) {
    start = start || 0;
    size = size || points.length;

    if ( size > 2 ) {
      var side1 = splitCurve(points, s, start, size-1 )[0];
      var side2 = splitCurve(points, s, start+1, size-1 )[1];
    } else {
      var side1 = [points[start]], side2 = [points[start+1]];
    }

    var p1 = side1[side1.length-1], p2 = side2[0];
    var midpoint = [vadd( p1, vscale( vsub(p2, p1), s) )];
    return [side1.concat(midpoint), midpoint.concat(side2)];
  }

  function Path( start ) {
    var segments = [], points = [start], path = {}, totalWeight = 0;

    path.curve = function(ps, weight) {
      if (weight == null) weight = 1;

      for (var i=0,point;point = ps[i];i++) points.push(point);

      segments.push({s:points.length-ps.length-1, o:ps.length+1, w:weight});

      totalWeight += weight;
      return path;
    }

    path.line = function(point, weight) {
      return path.curve([point], weight);
    }

    path.spline = function(ps, weight) {
      if (points.length < 2) path.curve(ps);

      var p1 = points[points.length-1];
      var p2 = points[points.length-2];
      var ctl = vadd( p1, vscale( vsub(p2,p1), -1) );

      return path.curve( [ctl].concat(ps), weight);
    }

    // Given a control point and an end point, add a cubic bezier curve that
    // approximates an arc of the circle that touches the current point and
    // the end point such that the tangents at both those points intersect at
    // the control point. This only works if the control point is equidistant
    // from the start and end points. If it isn't, a new control point that is
    // both equidistant and lies on the line between the start point and the
    // given control point will be calculated.
    path.arc = function(ctl, end, weight) {
      var pts = arcToCubic(points[points.length-1], ctl, end);
      return path.curve(pts.slice(1), weight);
    }

    path.vertices = function(numVertices, step) {
      return function(vertexSink) {
        var allDivisions = Math.max(numVertices, segments.length + 1);
        var remainingDivisions = allDivisions - 1 - segments.length;
        var remainingWeight = totalWeight;
        var index = 0;

        for (var i=0,s; s = segments[i]; i++) {
          var divisions = i >= segments.length - 1 ? remainingDivisions + 1 : (s.w * remainingDivisions / remainingWeight + 1)|0;
          for (var j=0; j<divisions; j++) {
            vertexSink( new Vertex(resolveCurve(points,j/divisions, s.s, s.o), step, (index++)/(allDivisions-1)) );
          }

          remainingDivisions -= divisions - 1;
          remainingWeight -= s.w;
        }

        vertexSink( new Vertex(points[points.length-1], step, 1) );
      }
    }

    path.startPoint = function() {
      return points[0];
    }

    // the given visitor function ([point], Float ->) will be passed an array of control points
    // and the weight for each segment in this path
    path.visitSegments = function(visitor) {
      for (var i=0, s; s = segments[i]; i++)
        visitor(points.slice(s.s+1,s.s+s.o), s.w)
    }

    // Create a new curve formed from the start point and all the segments of this path, and
    // the segments of the given path. The start point of the given path will be dropped.
    path.concat = function(secondPath) {
      var out = Path(points[0]);
      path.visitSegments(out.curve);
      secondPath.visitSegments(out.curve);
      return out;
    }

    path.transform = function(pointTransform) {
      var path = Path(pointTransform(points[0]));
      for (var i=0, s; s = segments[i]; i++)
        path = path.curve(points.slice(s.s+1,s.s+s.o).map(pointTransform),s.w);

      return path;
    }

    path.reverse = function() {
      var path = Path(points[points.length-1]);
      for (var i=segments.length-1, s; s = segments[i]; i--)
        path = path.curve(points.slice(s.s,s.s+s.o-1).reverse(),s.w);

      return path;
    }

    path.split = function(step) {
      var path1 = Path(points[0]);
      for (var i=0, s, weight=0; s = segments[i]; i++) {
        if ( (weight + s.w) / totalWeight >= step )
          break;
        path1 = path1.curve(points.slice(s.s+1,s.s+s.o),s.w);
        weight += s.w;
      }
      var split = splitCurve( points, (step*totalWeight-weight) / s.w, s.s, s.o );
      path1 = path1.curve( split[0].slice(1), step*totalWeight-weight );
      var path2 = Path(split[1][0]).curve(split[1].slice(1), s.w - step*totalWeight + weight);

      for (var j=i+1, s; s = segments[j]; j++)
        path2 = path2.curve(points.slice(s.s+1,s.s+s.o),s.w);

      return [path1,path2];
    }

    path.vertexAt = function(step) {
      for (var i=0, s, weight=0; s = segments[i]; i++) {
        if ( (weight + s.w) / totalWeight >= step )
          break;
        weight += s.w;
      }
      return new Vertex(resolveCurve( points, (step*totalWeight-weight) / s.w, s.s, s.o ), step, step );
    }

    path.stepSink = function(transformStep) {
      return function(vertexSink) {
        return function(step) {
          var remainingWeight = step;

          for (var i=0,s; s = segments[i]; i++) {
            var segmentWeight = s.w / totalWeight;
            if (remainingWeight - segmentWeight <= 0.00001) {
              return vertexSink( new Vertex(resolveCurve(points,remainingWeight, s.s, s.o), transformStep, step) );
            } else {
              remainingWeight -= segmentWeight;
            }
          }
        }
      }
    }

    // Given a predicate (Vertex -> Bool), find the first point on the path such that the predicate takes
    // different values on either side of it. The optional tolerance specifies how close (in parameter)
    // space the returned point must be to the border. Supplying a divisions parameter will start the
    // algorithm by dividing the curve that many times to avoid missing transitions (the default is 4).
    // The transform and ribstep of the returned vertex (and each vertex supplied to the predicate) will
    // be its parameter along the curve.
    path.vertexBordering = function(predicate, tolerance, initialDivisions) {
      tolerance = tolerance || .00001;
      divisions = initialDivisions || 4;
      var predicateValue = predicate(path.vertexAt(0));
      var step = 1 / (divisions-1);
      for (var parameter=0; parameter <= 1-step; parameter += step) {
        var testVertex = path.vertexAt(parameter+step);
        if (predicate(testVertex) != predicateValue) {
          while (step > tolerance) {
            step /= 2;
            testVertex = path.vertexAt(parameter+step);
            if (predicate(testVertex) == predicateValue)
              parameter += step;
          }
          return testVertex;
        }
      }
    }

    return path;
  }

  // STEP
  function step(iterations) {
    return function(stepSink) {
      for (var i=0; i < iterations; i++) stepSink( i / (iterations - 1) );
    }
  }

  // LIFT

  // Manage a series of geometry computations.
  // lift should be called as
  //   lift(source).generate(generator1,generator2,...).render(renderer)
  // Where source is an ASink ->, a generator1 is BSink -> ASink, generator2 is CSink -> BSink, etc.
  // and renderer is a FaceSink. The last generator argument must accept a FaceSink.
  // Nothing will be executed until render is called.
  function lift(source) { return {
    generate: function() {
      var generators = Array.prototype.slice.call(arguments,0); return {
      render: function(renderer) {
        source(generators.reverse().reduce(function(m,generator) { return generator(m); }, renderer));
      } }
    } }
  }

  // VERTEX GENERATORS

  // Given a list of 3d points, generated a vertex for each point per step
  // (points need to be transformed)
  function vertices(points) {
    return function(vertexSink) {
      return function(step) {
        for (var i=0,p; p=points[i]; i++)
          vertexSink( new Vertex(p, step, i / (points.length-1)) );
      }
    }
  }

  // Convert steps into vertices using a given function.
  function parametric(f, stepper) {
    return function(vertexSink) {
      return function( transformStep ) {
        stepper( function(ribStep) {
          vertexSink( new Vertex(f(ribStep, transformStep), transformStep, ribStep) );
        })
      }
    }
  }

  // Simplifies step generators
  function vertexGenerator( f ) {
    return function(vertexSink) {
      return function(step) {
        f(step,vertexSink)
      }
    }
  }

  // Convert a path into a manifold using a path generator function that converts vertices
  // on the path into new paths.
  function PathParameterized(path, transformSteps, ribSteps) {
    return function(generator) {
      return function(vertexSink) {
        path.vertices(transformSteps, 0)(function(vertex) {
          generator(vertex).vertices(ribSteps,vertex.ribStep)(vertexSink);
        });
      }
    }
  }

  // Given a pathGenerator (Vertex -> Path), creates an Vertex -> VertexSink stream
  // that generates ribSteps vertices along a each path generated from each
  // input vertex.
  // The transform step of each output vertex will be the rib step of each input vertex.
  function VertexParameterizedPath(pathGenerator, ribSteps) {
    return function(vertexSink) {
      return function(vertex) {
        pathGenerator(vertex).vertices(ribSteps,vertex.ribStep)(vertexSink);
      }
    }
  }

  // creates a VertexSource that, when Given a list of paths and a number of steps,
  // emits a point on each path.
  function PathSource(paths, steps) {
    return function(vertexSink) {
      for (var i=0, step=0; i < steps; i++, step=i/(steps-1)) {
        for (var j=0, path; path=paths[j]; j++) {
          var vertex = path.vertexAt( step );
          vertex.transformStep = step;
          vertex.ribStep = j / (paths.length - 1);
          vertexSink( vertex );
        }
      }
    }
  }

  // Given a list of paths, emit n lists of vertices along those paths, where n is
  // the given number of steps ([Path] -> [Vertex]). The points on each path will
  // be linearly spaced in parameter space (meaning curves will NOT be linearly spaced
  // in 3d space). Both the transform and rib steps of each vertex will be the step
  // parameter used to generate the vertex.
  function RibTransform(paths, steps) {
    return function(verticesSink) {
      for (var i=0; i<steps; i++)
        verticesSink( paths.map(function(path) { return path.vertexAt( i/(steps-1) ); }) )
    }
  }

  // Given a rib transform ([Vertex] -> A) and a multi-vertex path Generator ([Vertex] -> Path),
  // emit vertices for a manifold.
  function RibParameterizedPath( ribTransform, pathGenerator, ribSteps ) {
    return function(vertexSink) {
      ribTransform( function(vertices) {
        pathGenerator(vertices).vertices(ribSteps,vertices[0].ribStep)(vertexSink);
      })
    }
  }

  // A vertex transformer that sequences a list of vertex transformers.
  // Sequences look like:
  // [[3, transformerA], [5, transformerB], ...]
  // This will send the first 3 vertices from the sequencer source to transformerA, the
  // next 5 vertices to transformer, etc. The ribStep of each input vertex will be put
  // on the range [0 1] so that each transformer receives a full range over all of its
  // steps. Each vertex emitted to the sequence's sink will have it's transform step
  // translated using the opposite transform.
  function Sequencer(sequence) {
    return function(vertexSink) {

      function mapTransformStep(vertex) {
        vertex.transformStep = Math.min(1, step / (totalSteps-1) );
        vertexSink( vertex );
      }

      sequence = sequence.map(function(s) { return [s[0], s[1]( mapTransformStep )]; });

      var current = sequence[0];
      var rest = sequence.slice(1);
      var step = 0, sequenceStep = 0;
      var totalSteps = sequence.reduce(function(m,s) {return s[0]+m;},0);

      return function(vertex) {

        while (current && step - sequenceStep >= current[0]) {
          sequenceStep += current[0];
          current = rest.splice(0,1)[0];
        }

        if (current) {
          var ribStep = current[0] == 1 ? 1 : (step - sequenceStep) / (current[0]-1);
          current[1](new Vertex(vertex, ribStep, vertex.ribStep));
          step++;
        }
      }
    }
  }

  // Takes a list of weights and vertexSource and combines them into a single source.
  // [[Float, vertexSink ->]] -> vertexSink
  // The transform step of each source will be altered according to the weight given to
  // each source. The sources are assumed to run synchronously.
  // The spacing of the transform step is an optional second parameter.
  function Stage(stages, spacing) {
    var totalWeight = stages.reduce(function(m,stage) { return m+stage[0]; }, 0);
    var spacing = spacing || totalWeight * .01;
    totalWeight += spacing * (stages.length - 1);

    return function(vertexSink) {
      var cumulativeWeight = 0;
      stages.forEach(function(stage) {
        stage[1](function(vertex) {
          vertex.transformStep = Math.min(1, (vertex.transformStep * stage[0] + cumulativeWeight) / totalWeight );
          vertexSink(vertex);
        })
        cumulativeWeight += stage[0] + spacing;
      });
    }
  }

  // Create a vertext generator that converts a single vertex into a circle of vertices
  // such that the circle passes through the vertex and its center lies on the center normal.
  function CircleRib( steps, centerNormal, centerOffset, phase ) {
    centerNormal = vnorm(centerNormal);
    centerOffset = centerOffset || [0,0,0];
    centerOffset = vsub(centerOffset, vscale(centerNormal, vdot(centerOffset,centerNormal)));
    phase = phase || 0;
    return function(vertexSink) {
      return function(vertex) {
        var center = vadd(vscale(centerNormal, vdot(centerNormal,vertex)), centerOffset);
        var a = vsub(vertex,center);
        var radius = vlength(a);
        if (radius === 0)
          return vertexSink( new Vertex(center, vertex.ribStep, 1) );

        var b = vcross(centerNormal,a);
        for (var i=0; i<steps; i++) {
          var point = vadd(vadd(vscale(a,Math.cos(phase+2*i*Math.PI/(steps))),vscale(b,-Math.sin(phase+2*i*Math.PI/(steps)))),center);
          vertexSink( new Vertex(point, vertex.ribStep, i/(steps-1)) );
        }
      }
    }
  }

  // TRANSFORM
  function translate(translations) {
    return function( vertexSink ) {
      return function( vertex ) {
        var tIndex = parseInt(vertex.transformStep * (translations.length-1))
        var translate = vertex.transformStep==1 ?
          translations[tIndex] :
          vinterp( translations[tIndex], translations[tIndex+1], vertex.transformStep-tIndex)
        vertexSink( new Vertex(vadd(translate,vertex), vertex.transformStep, vertex.ribStep ) );
      }
    }
  }

  // TESSELATE
  function skin( faceSink ) {
    var lastRib;
    var nextRib = [];
    var bIndex = 0;
    var lastTransformStep = 0;
    return function skinVertexSink( vertex ) {
      if ( vertex.transformStep > lastTransformStep ) {
        lastRib = nextRib;
        nextRib = [];
        bIndex=0;
        lastTransformStep = vertex.transformStep;
      }

      if (lastRib) {
        var blVertex = lastRib[bIndex];
        var brVertex = lastRib[bIndex+1];

        if (nextRib.length) {
          faceSink( [blVertex, nextRib[nextRib.length-1], vertex] );
        }

        while ( brVertex &&  vertex.ribStep > blVertex.ribStep + (brVertex.ribStep-blVertex.ribStep)/2 ) {
          faceSink( [blVertex, vertex, brVertex] );
          blVertex = brVertex;
          brVertex = lastRib[++bIndex + 1];
        }
      }

      nextRib.push(vertex)
    }
  }

  function skinRibs( rib1, rib2, faceSink ) {
    var bIndex = 0;
    rib2.forEach(function(vertex,i) {
      var blVertex = rib1[bIndex];
      var brVertex = rib1[bIndex+1];

      if (i>0 && rib2.length > i) {
        console.log(printVertex(blVertex), printVertex(rib2[i-1]), printVertex(vertex))
        faceSink( [blVertex, rib2[i-1], vertex] );
      }

      while ( brVertex &&  vertex.ribStep > blVertex.ribStep + (brVertex.ribStep-blVertex.ribStep)/2 ) {
        faceSink( [blVertex, vertex, brVertex] );
        console.log('while', printVertex(blVertex), printVertex(vertex), printVertex(brVertex))
        blVertex = brVertex;
        brVertex = rib1[++bIndex + 1];
      }
    })
  }

  function facers(facer1, facer2, etc) {
    var facers = Array.prototype.slice.call(arguments,0);
    return function facersFacer( faceSink ) {
      var vertexSinks = facers.map( function(f) { return f(faceSink); } );
      return function facersVertexSink( vertex ) {
        vertexSinks.forEach(function(vs) {vs(vertex)});
      }
    }
  }

  function closeEdge( faceSink ) {
    var bottomFirstInRib, bottomLastInRib;
    var topFirstInRib, topLastInRib;
    return function wrapEdgeVertexSink( vertex ) {
      if ( ! topFirstInRib )
        topFirstInRib = vertex;
      else if ( vertex.transformStep == topFirstInRib.transformStep )
        topLastInRib = vertex;

      if ( vertex.ribStep == 1 || vertex.transformStep > topFirstInRib.transformStep ) {
        if ( bottomFirstInRib ) {
          if ( topLastInRib && bottomLastInRib )
            faceSink( [bottomLastInRib, topLastInRib, topFirstInRib] );

          if ( bottomLastInRib )
            faceSink( [bottomLastInRib, topFirstInRib, bottomFirstInRib] );
          else if ( topLastInRib )
            faceSink( [bottomFirstInRib, topLastInRib, topFirstInRib] );
        }

        bottomFirstInRib = topFirstInRib;
        bottomLastInRib = topLastInRib;
        topFirstInRib = vertex.transformStep > topFirstInRib.transformStep ? vertex : null;
        topLastInRib = null;
      }
    }
  }

  function capBottom( faceSink ) {
    var rib = [];
    return function capBottomVertexSink( vertex ) {
      if ( vertex.transformStep === 0 )
        return rib.push(vertex);

      if (rib) {
        tesselate(rib, reverseFaceSink(faceSink) );
        rib = null;
      }
    }
  }

  function capTop( faceSink ) {
    var rib = [];
    return function capTopVertexSink( vertex ) {
      if ( vertex.transformStep < 1 ) return;

      rib.push(vertex);

      if (vertex.ribStep == 1)
        tesselate(rib, faceSink );
    }
  }

  function capTube( faceSink ) {
    var outerRib = [];
    var innerRib = []
    return function capTubeVertexSink( vertex ) {
      if ( vertex.transformStep == 0 ) {
        outerRib.push(vertex);
      } else if ( vertex.transformStep == 1 ) {
        innerRib.push(vertex);
        if (vertex.ribStep == 1) {
          // add outer rib point again to close path
          outerRib.push(outerRib[0]);

          // start inner rib at closest point
          // var start = outerRib[0];
          // var closest = innerRib.reduce(function(m,p,i) { return vdist(start,innerRib[m]) > vdist(start,p) ? i-1 : m; }, 0);
          // innerRib = innerRib.slice(closest).concat( innerRib.slice(0,closest) );

          // // add inner rib point again to close path
          // innerRib.push(innerRib[0]);
          // innerRib.reverse();
          // tesselate( outerRib.concat(innerRib), reverseFaceSink(faceSink));
          skinRibs( innerRib, outerRib, faceSink );
        }
      }
    }
  }

  function capTubeBottom() {
    var rib = [];
    var op = outerManifoldStart;

    return function( faceSink ) {
      function outerManifoldStart(vertex) {
        if ( vertex.transformStep > 0 ) return outerManifoldMid;
        rib.push(vertex);
        return outerManifoldStart;
      }

      function outerManifoldMid(vertex) {
        return vertex.transformStep < 1 && vertex.ribStep < 1 ? outerManifoldMid : innerManifoldStart;
      }

      function innerManifoldStart(vertex) {
        if ( vertex.transformStep > 0 ) {
          tesselate(rib, faceSink );
          return noop;
        }
        rib.push(vertex);
        return outerManifoldStart;
      }

      function noop() { }

      return function( vertex ) { op = op(vertex); }
    }
  }

  function capTubeTop() {
    var rib = [];
    var op = outerManifoldStart;

    return function( faceSink ) {
      function outerManifoldStart(vertex) {
        if ( vertex.transformStep < 1 ) return outerManifoldStart;
        rib.push(vertex);
        return outerManifoldEnd;
      }

      function outerManifoldEnd(vertex) {
        rib.push(vertex);
        if ( vertex.ribStep < 1 ) return outerManifoldEnd;
        return innerManifoldStart;
      }

      function innerManifoldStart(vertex) {
        rib.push(vertex);
        if ( vertex.ribStep < 1 ) return innerManifoldStart;

        tesselate(rib, faceSink );
        return noop;
      }

      function noop() { }

      return function( vertex ) { op = op(vertex); }
    }
  }

  function debugFacer( faceSink ) {
    return function( vertex ) { console.log(vertexString(vertex), vertex.transformStep, vertex.ribStep); };
  }

  // FACER TRANSFORMS

  // Facer -> Facer
  function reverse(facer) {
    return function( faceSink ) {
      return facer( reverseFaceSink(faceSink) )
    }
  }

  // FaceSink -> FaceSink
  function reverseFaceSink( faceSink ) {
    return function(face) {
      faceSink( [face[2],face[1],face[0]] )
    }
  }

  function vertexString(v) {return v?'{'+v.x+','+v.y+','+v.z+'} ':'null'}
  function faceString(face) { return face.map(vertexString); }

  // RENDER
  function ThreeJSRenderer() {
    var geometry = new THREE.Geometry();
    geometry._manifoldIdMap = {};
    return {
      renderer : function( face ) {
        var index0 = saveThreeJSVertex( geometry, face[0] );
        var index1 = saveThreeJSVertex( geometry, face[1] );
        var index2 = saveThreeJSVertex( geometry, face[2] );
        geometry.faces.push( new THREE.Face3(index0, index1, index2) );
      },
      geometry: geometry
    }
  }

  function CSGRenderer() {
    var polygons = [];
    return {
      renderer : function( face ) {
        var faceNormal = vnorm(vcross(vsub(face[0],face[1]),vsub(face[2],face[1])));
        var vertices = [
          new CSG.Vertex( face[0], faceNormal ),
          new CSG.Vertex( face[1], faceNormal ),
          new CSG.Vertex( face[2], faceNormal )];
        polygons.push( new CSG.Polygon( vertices ) );
      },
      csgObject: function() { return CSG.fromPolygons(polygons); }
    }
  }

  function saveThreeJSVertex(geometry, vertex) {
    var location = geometry._manifoldIdMap[vertex.id];
    if ( location == null ) {
      location = geometry._manifoldIdMap[vertex.id] = geometry.vertices.length;
      geometry.vertices.push( new THREE.Vector3( vertex[0], vertex[1], vertex[2] ) );
    }
    return location;
  }

  function STLRenderer(){
    var doc = 'solid pixel';
    return {
      renderer : function( face ) {
        var normal = vnorm( vcross( vsub(face[1],face[0]), vsub(face[2],face[0])) );
        doc += "facet normal " + normal.x + " " + normal.y + " " + normal.z + " \n";
        doc += "outer loop \n";
        doc += "vertex " + face[0].x + " " + face[0].y + " " + face[0].z + " \n";
        doc += "vertex " + face[1].x + " " + face[1].y + " " + face[1].z + " \n";
        doc += "vertex " + face[2].x + " " + face[2].y + " " + face[2].z + " \n";
        doc += "endloop \n";
        doc += "endfacet \n";
      },
      doc : function() { return doc + 'endsolid' }
    }
  }


  // MATH

  function vadd(a,v) { return [a[0]+v[0], a[1]+v[1], a[2]+v[2]]; }
  function vsub(a,v) { return [a[0]-v[0], a[1]-v[1], a[2]-v[2]]; }
  function vscale(a,c) { return [a[0]*c, a[1]*c, a[2]*c]; }
  function vdot(a,v) { return a[0]*v[0] + a[1]*v[1] + a[2]*v[2]; }
  function vcross(a,v) { return [a[1]*v[2] - a[2]*v[1], a[2]*v[0] - a[0]*v[2], a[0]*v[1] - a[1]*v[0]]; }
  function vlength(v) { return Math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }
  function vdist(a, v) { return vlength(vsub(a,v)); }
  function vnorm(v) { var l=vlength(v); return l > 0 ? [v[0]/l, v[1]/l, v[2]/l] : v; }
  function vinterp(a,b,c) { return vadd(a,vscale(vsub(b,a),c)); }

  function vectorAverage(vs) { return vs.length ? vscale( vs.reduce(vadd,[0,0,0]), 1/vs.length ) : [0,0,0]; }

  var Q = {
    quaternion : function(vec) { return [0, vec[0], vec[1], vec[2]]; },
    rotation : function( a, vec ) { var a2=a/2, sina2=Math.sin(a2);
               return [Math.cos(a2), vec[0]*sina2, vec[1]*sina2, vec[2]*sina2 ]; },
    safeRotation : function( a, vec ) { return Q.rotation(a,vnorm( vec )); },
    rotate : function( vec, rot ) { return Q.mul(Q.mul(rot,Q.quaternion(vec)),Q.conjugate(rot)).slice(1); },
    conjugate : function(q) { return [q[0],-q[1],-q[2],-q[3]] },
    mul : function(q1,q2) { return [
          q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3],
          q1[0]*q2[1] + q1[1]*q2[0] - q1[2]*q2[3] + q1[3]*q2[2],
          q1[0]*q2[2] + q1[1]*q2[3] + q1[2]*q2[0] - q1[3]*q2[1],
          q1[0]*q2[3] - q1[1]*q2[2] + q1[2]*q2[1] + q1[3]*q2[0] ]; },
  }

  function loopMeanNormal(vs) {
    function v(i) { return vs[i%vs.length] }
    return vectorAverage(vs.map(function(v0,i){return vcross(vsub(v0,v(i+1)), vsub(v(i+2), v(i+1)))}))
  }

  function intersect(v1, v2, v3, v4) {
    var det = (v1.x-v2.x)*(v3.y-v4.y) - (v1.y-v2.y)*(v3.x-v4.x);
    return !det ? null : {
      x: ((v1.x*v2.y-v2.x*v1.y)*(v3.x-v4.x) - (v3.x*v4.y-v4.x*v3.y)*(v1.x-v2.x))/det,
      y: ((v1.x*v2.y-v2.x*v1.y)*(v3.y-v4.y) - (v3.x*v4.y-v4.x*v3.y)*(v1.y-v2.y))/det
    };
  }

  function doesIntersect(v1, v2, v3, v4) {
    var p = intersect(v1,v2,v3,v4);
    if (!p) return false;
    return p.x >= Math.min(v1.x,v2.x) && p.x >= Math.min(v3.x,v4.x)
      && p.x <= Math.max(v1.x,v2.x) && p.x <= Math.max(v3.x,v4.x)
      && p.y >= Math.min(v1.y,v2.y) && p.y >= Math.min(v3.y,v4.y)
      && p.y <= Math.max(v1.y,v2.y) && p.y <= Math.max(v3.y,v4.y);
  }

  function isSelfIntersecting( p1, p2, vertices ) {
    for (var i=0,v1,v2; v1=vertices[i], v2=vertices[i+1]; i++)
      if ( doesIntersect(p1,p2,v1,v2) ) return true;
    return false;
  }

  function projectIntoPlane( p, planePoint, planeNormal ) {
    // p - pp - (((p-pp) . pn) * pn)
    var pTranslated = vsub( p, planePoint );
    return vsub( pTranslated, vscale( planeNormal, vdot( pTranslated, planeNormal )))
  }

  function loop( vertices, index ) {
    return vertices[ (index + vertices.length) % vertices.length ];
  }

  function arcToCubic(start, ctl, end) {
    // find the control point equidstant between start and the end point.
    var s = (2*vdot(start,end) - vdot(start,start) - vdot(end,end)) / (2*vdot(vsub(ctl,start), vsub(start,end)));
    ctl = vadd(start, vscale(vsub(ctl,start), s) );

    // cross the tangent lines to get a normal to the plane of the circle.
    var planeNormal = vnorm(vcross(vsub(end,ctl), vsub(start,ctl)));

    // cross each tangent with the normal to find the direction vectors from the points to the center.
    var d1 = vnorm( vcross( planeNormal, vsub(ctl,start) ) );
    var d2 = vnorm( vcross( planeNormal, vsub(end,ctl) ) );

    // The intersection of the direction vectors gives the center.
    s = vlength(vsub(end,start)) / vlength(vsub(d1,d2));
    var center = vadd(start, vscale(d1,s));

    // Cross the direction vectors to get the angle betwen the points
    var angle = Math.asin(Math.min(1,vlength(vcross(d1,d2))));

    // Compute length of ctl points (derived from http://itc.ktu.lt/itc354/Riskus354.pdf)
    var cosa = Math.cos(angle/2);
    var sina = Math.sin(angle/2);
    var radius = vdist(start, center);
    var ctlLength = radius * vlength([(4-cosa)/3 - cosa, (1-cosa)*(3-cosa)/(3*sina) - sina,0]);

    // if angle > pi/2, invert ctlLength so control points keep getting longer
    if (vdot(d1,d2) < 0)
      ctlLength = radius*(8*(Math.sqrt(2)-1)/3) - ctlLength;

    // control points on tangents, ctlLength away from start and end points
    var ctl1 = vadd( start, vscale( vnorm(vsub(ctl,start)) ,ctlLength ) );
    var ctl2 = vadd( end, vscale( vnorm(vsub(ctl,end)) ,ctlLength ) );

    return [start, ctl1, ctl2, end];
  }

  // TESSLATION

  // Assume 2d with l2.x >= v.x >= l1.x
  function vertexIsAboveLine( v, l1, l2 ) {
    return (v[0] == l1[0] && v[1] > l1[1] && (l2[0] > l1[0] || v[1] > l2[1]))
          || (v[1]-l1[1])/(v[0]-l1[0]) > (l2[1]-l1[1])/(l2[0]-l1[0]);
  }

  function isConvex(v1, v2, v3) {
    return vcross(vsub(v1,v2),vsub(v3,v2))[2] > 0;
  }

  function printVertex(v) { return '('+v.id+':'+v.x.toFixed(2)+','+v.y.toFixed(2)+','+v.z.toFixed(2)+')'}

  function MonotonePolygon( first, lowerVertices ) {
    var upper = first ? [first] : [];
    var lower = lowerVertices || [];
    var mergePolygon;

    function addUpper( v, faceSink ) {
      while ( upper.length > 1 && isConvex(upper[upper.length-2],upper[upper.length-1],v) ) {
        faceSink( [upper[upper.length-1], v, upper[upper.length-2]] );
        upper.pop();
      }

      if (upper.length == 1 && lower.length && isConvex(lower[0],upper[0],v))
        faceSink( [lower[0], upper.pop(), v] );

      for (var i=0, vl1, vl2; vl1 = lower[i], vl2=lower[i+1]; i++)
        faceSink( [vl1, v, vl2] );

      upper.push(v);
      lower = lower.slice(lower.length-1,lower.length);
    }

    function addLower( v, faceSink ) {
      while ( lower.length > 1 && isConvex(v,lower[lower.length-1],lower[lower.length-2]) ) {
        faceSink( [v, lower[lower.length-1], lower[lower.length-2]] );
        lower.pop();
      }

      if (upper.length && lower.length == 1 && isConvex(v,lower[0],upper[0]))
        faceSink( [lower.pop(), upper[0], v] );

      for (var i=0, vl1, vl2; vl1 = upper[i], vl2=upper[i+1]; i++)
        faceSink( [vl2, v, vl1] );

      lower.push(v);
      upper = upper.slice(upper.length-1,upper.length);
    }


    function lastLower() {
      return lower[lower.length-1] || upper[0]
    }

    function aboveBottom(v, vertices) {
      var vl1 = lastLower();
      var vl2 = loop( vertices, vl1.id - 1 );
      return vertexIsAboveLine(v, vl1, vl2 )
    }

    function isBottom(v, vertices) {
      return v === loop( vertices, lastLower().id - 1 );
    }

    function getUpper() { return upper; }
    function getLower() { return lower; }

    function attemptAdd( v, vertices, faceSink ) {
      if (debug)
        console.log('vertex',v, upper.map(printVertex).join(','), lower.map(printVertex).join(','));
      var vu1 = upper[upper.length-1];
      if ( ! vu1 ) {
        addUpper( v, faceSink );
        return "TOP";
      }

      var vu2 = loop( vertices, vu1.id+1 );
      if ( v === vu2 ) {
        addUpper( v, faceSink );
        if (mergePolygon) {
          mergePolygon.addUpper(v, faceSink);
          upper = mergePolygon.getUpper();
          lower = mergePolygon.getLower();
          mergePolygon = null;
        }
        return isBottom(v,vertices) ? "DONE" : "TOP";
      }

      if ( mergePolygon ? mergePolygon.isBottom(v,vertices) : isBottom(v,vertices) ) {
        addLower( v, faceSink );
        if (mergePolygon) {
          mergePolygon.addLower( v, faceSink );
          mergePolygon = null;
        }
        var vl1 = lastLower();
        var vl2 = loop( vertices, vl1.id - 1 );
        return vl2.x >= vl1.x ? "BOTTOM" : "MERGE";
      }

      if ( vertexIsAboveLine(v, vu1, vu2 ) )
        return "ABOVE";
      else if ( mergePolygon ? mergePolygon.aboveBottom(v,vertices) : aboveBottom(v,vertices) )
        return "INSIDE";
      else
        return "BELOW";
    }

    function merge( polygon ) {
      mergePolygon = polygon;
    }

    function split( v, faceSink ) {
      var other = mergePolygon || MonotonePolygon( null, lower.length ? [lastLower()] : [upper[0]] );
      mergePolygon = null;

      if (debug) {
        console.log('SPLIT', upper.map(printVertex).join(','), lower.map(printVertex).join(','));
        console.log('OTHER', other.getUpper().map(printVertex).join(','), other.getLower().map(printVertex).join(','))
      }
      other.addUpper( v, faceSink );
      addLower( v, faceSink );

      return other;
    }

    return { addUpper:addUpper, addLower:addLower, attemptAdd:attemptAdd, merge:merge, split:split,
             isBottom:isBottom, aboveBottom:aboveBottom, getUpper:getUpper, getLower:getLower };
  }


  function tesselate( points3d, faceSink ) {
    // determine approximate plane through loop and map points onto it.
    var planeNormal = vnorm(loopMeanNormal(points3d));
    var planePoint = vectorAverage(points3d);

    var offNormal = vadd( planePoint, Math.abs(planeNormal[0] < .5) ? [1,0,0] : [0,1,0] );
    var axis1 = vnorm(projectIntoPlane(offNormal ,planePoint, planeNormal));
    var axis2 = vcross( planeNormal, axis1 );

    var vertices2d = points3d.map( function(p,index) {
      var pointInPlane = projectIntoPlane(p, planePoint, planeNormal);
      return new Vertex([vdot(pointInPlane,axis1), vdot(pointInPlane,axis2),0],0,0,index);
    })

    // tesselate 2d loop
    tesselate2d( vertices2d, function(face) {
      faceSink([ points3d[face[2].id], points3d[face[1].id], points3d[face[0].id] ]);
    })
  }

  function tesselate2d( vertices, faceSink ) {
    var sortedVertices = vertices.slice().sort(function(a,b) { return a.x - b.x || b.y - a.y; });
    var monotones = [MonotonePolygon(sortedVertices[0])];

    for (var i=1, v; v = sortedVertices[i]; i++) {
      var result;
      for (var j=0,mp; mp = monotones[j]; j++) {
        result = mp.attemptAdd(v, vertices, faceSink);
        if (result !== 'BELOW') break;
      }
      if (debug)
        console.log( "TESSELATE", j, result )

      if ( result === 'ABOVE' )
        monotones.splice( j, 0, MonotonePolygon(v) );
      else if ( result === 'INSIDE' )
        monotones.splice( j+1, 0, mp.split(v, faceSink) );
      else if ( result === 'BELOW' )
        monotones.push( MonotonePolygon(v) );
      else if ( result === 'MERGE' ) {
        monotones[j+1].attemptAdd( v, vertices, faceSink );
        monotones[j].merge(monotones[j+1])
        monotones.splice(j+1,1)
      } else if (result === 'DONE' )
        monotones.splice(j, 1);
    }
  }

  var all = {
      vadd:vadd, vsub:vsub, vscale:vscale, vdot:vdot, vcross: vcross, vlength:vlength, vnorm:vnorm,
      Q:Q, step:step,
      Path:Path, PathParameterized:PathParameterized, VertexParameterizedPath:VertexParameterizedPath,
      RibTransform:RibTransform, RibParameterizedPath:RibParameterizedPath, PathSource:PathSource,
      lift:lift, CircleRib:CircleRib, Sequencer:Sequencer, Stage:Stage,
      Vertex:Vertex, vertices:vertices, parametric:parametric, vertexGenerator:vertexGenerator,
      translate:translate,MonotonePolygon:MonotonePolygon, tesselate2d:tesselate2d,
      skin:skin, facers:facers, closeEdge:closeEdge, capBottom:capBottom, capTop:capTop, capTube:capTube,
      debugFacer:debugFacer, reverse:reverse, arcToCubic:arcToCubic,
      ThreeJSRenderer:ThreeJSRenderer, CSGRenderer:CSGRenderer, STLRenderer:STLRenderer
  };
  for (var k in all) context[k] = all[k];

})(typeof window != 'undefined' ? (window.manifold={}) : (typeof self != 'undefined' ? (self.manifold={}) : exports ));
