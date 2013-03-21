"use strict"

var numeric = require("numeric")
var top = require("simplicial-complex")
var differential = require("differential")
var cooriented = require("cooriented")

var EPSILON = 1e-10

//Higher dimensional network flow:
//
//  Given a simplicial complex, let d be the differential, e and f be chains and alpha/w_i be cochains
//
//Then the n-dimensional flow problem is the following linear program:
//
//  Minimize alpha(f)
//  s.t.
//      d(f) = e
//      0 <= w_i(f) <= c_i
//
function ndflow(cells, capacities, e_cells, e_weights, alpha) {

  //Want to put problem into form:
  //
  //  Minimize c . x
  //  s.t.
  //        Ax <= b
  //
  var c = alpha
  var A = []
  var b = []
  
  //Assemble conservation constraint:
  //
  //    d(f) = e
  //
  var D = differential(cells)
  var Dd = D.toDense()
  var ep = numeric.rep([D.boundaryCells.length], 0.0)
  for(var i=0; i<e_cells.length; ++i) {
    var idx = top.findCell(D.boundaryCells, e_cells[i])
    if(idx < 0) {
      return {
        cells: [],
        weights: []
      };
    }
    var orientation = cooriented(D.boundaryCells[idx], e_cells[i])
    var v = -e_weights[i] * orientation
    if(v < 0) {
      ep[idx] = -v
    } else {
      Dd[idx] = numeric.neg(Dd[idx])
      ep[idx] = v
    }
  }
  A = A.concat(Dd).concat(numeric.neg(Dd))
  b = b.concat(ep).concat(numeric.add(EPSILON, numeric.neg(ep)))
  
  //Assemble capacity constraints:
  //
  //    0 <= w_i(f) <= c_i
  //
  for(var i=0; i<cells.length; ++i) {
    //Add 0 <= w_i(f) constraint
    var w_i = numeric.rep([cells.length], 0.0)
    w_i[i] = -1.0
    A.push(w_i)
    b.push(capacities[i])
    
    //Add w_i(f) <= c_i(f)
    w_i = numeric.rep([cells.length], 0.0)
    w_i[i] = 1.0
    A.push(w_i)
    b.push(capacities[i])
  }

  //Solve!
  var result = numeric.solveLP(c, A, b)
  
  //Unpack into a chain
  var r_cells = []
  var r_weights = []
  if(result.solution) {
    for(var i=0; i<cells.length; ++i) {
      if(Math.abs(result.solution[i]) > EPSILON) {
        r_cells.push(cells[i])
        r_weights.push(result.solution[i])
      }
    }
  }
  return {
    cells: r_cells,
    weights: r_weights
  }
}
module.exports = ndflow
