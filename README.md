ndflow
======
Network flow in the nth-dimension!

Usage
=====
First install like this:

    npm install ndflow
    
Then use it as follows:

```javascript
require("ndflow")(
    [[0,1,2], [1,2,3]],           //Simplicial complex
    [1, 1],                       //Capacities
    [[0,1], [1,3], [3,2], [2,0]], //Boundary circulation
    [1,1,1,1],                    //Boundary weights
    [1,1])                        //Cost function
```

### `require("ndflow")(cells, capacities, boundary_cells, boundary_weights, cost)`
Solves for the flow from the boundary cells that minimizes the cost.  This is the same thing as solving the linear program:


```
Minimize      cost(flow)
Subject to    D(flow) = boundary_cells * boundary_weights
              0 <= flow <= capacities
```

In the case where the cell complex is one dimensional, this is the same as solving minimum cost network flow for a given capacity.  Binary searching on the scale of the capacities lets one compute a maximum/min-cost flow.

Credits
=======
(c) 2013 Mikola Lysenko. BSD License
