# -*- coding: utf-8 -*-
"""
*************** PATCH TEST *********************************


. Mesh olarak 'cube_in_cube.mesh' kullanıldı. Connectivity listte'ki en sondaki rakamlar elemanin group'unu gosterir.
. Right and Left Regions x<-5.9 ve x>5.9 olarak ayarlandı.
. x yönünde Predisplacement vermek için get_boundary function eklendi. 
. activeFibers çıkarıldı. Bu yüzden vf_matrix = 1.0 yapıldı.
. Right region boundary condition function ile uygulandı.



"""

import numpy as nm

from sfepy import data_dir
filename_mesh = data_dir + '/meshes/3d/cube_in_cube.mesh'
bc_value = 1.0
vf_matrix = 1.0
#vf_fibres1 = 0.8

#######################################################################################################
"""
.  'nls' : 'newton'  ->  string : nonlinear solver name 
.  'ls' : 'ls'       ->  string : linear solver name
.  'ts' : 'ts'       ->  string : time stepping solver name
.  'save_steps' : -1 ->  string : int, number of time steps when results should be saved and -1 for all time steps
.  'post_process_hook' : 'stress_strain'  -> string : a function to be called after each time step, 
                                                      used to update the results to be saved
                                                      
"""

options = {
    'nls' : 'newton',
    'ls' : 'ls',
    'ts' : 'ts',
    'save_steps' : -1,
    'post_process_hook' : 'stress_strain',
    'output_dir' : '/home/gkdnz/Desktop/PlotOutput/vtk',
}
#######################################################################################################
"""
.  displacement : field name
.  nm.float64   : datatype
.  3            : number of DOFs per node
.  Omega        : the name of region where the field is defined
.  1            : FE approximation order
"""

fields = {
    'displacement': (nm.float64, 3, 'Omega', 1),
}

#######################################################################################################
materials = {
    'solid' : ({
    
        'K'  : vf_matrix * 1e3, # bulk modulus, solid.K yazılarak çağrılır.
        'mu' : vf_matrix * 20e0, # shear modulus of neoHookean term, solid.mu yazılarak çağrılır.
    },),
    #'f1' : 'get_pars_fibres1', # equations kısmında get_pars_fibres1 yerine f1 kullanılacak.
}

#######################################################################################################
"""
def get_pars_fibres(ts, coors, mode=None, which=0, vf=1.0, **kwargs):
    
    Parameters
    ----------
    ts : TimeStepper
        Time stepping info.
    coors : array_like
        The physical domain coordinates where the parameters shound be defined.
    mode : 'qp' or 'special'
        Call mode.
    which : int
        Fibre system id.
    vf : float
        Fibre system volume fraction.
    
    
    
    Açıklamalar which = 0 için yazılmıştır.    
        
    fdir.shape =(3,1) --> fdir = array([[ 1.],
                                       [ 0.],
                                       [ 0.]])


    if mode != 'qp': return

    fmax = 1
    eps_opt = 0.01
    s = 0.2
    
    tt = ts.nt * 2.0 * nm.pi
    print "ts.nt: ", ts.nt;
    print "self.time: ",ts.time;
    print "t0: ",ts.t0
    print "t1: ",ts.t1
    print "dt: ",ts.dt
    print "step: ",ts.step
    if which == 0: # system 1
        fdir = nm.array([1.0, 0.0, 0.0], dtype=nm.float64) # fdir = array([1., 0., 0.])
        act = Factivation(ts.time)        
        #act = 0.5 * (1.0 + nm.sin(tt - (0.5 * nm.pi))) # pyfem'deki aktivaston function içindeki fa ile aynı.

    else:
        raise ValueError('unknown fibre system! (%d)' % which)

    fdir.shape = (3, 1) # coloumn matrix çevirir.(vector)
    fdir /= nm.linalg.norm(fdir) # nm.linalg.norm(fdir) = 1

    print act

    shape = (coors.shape[0], 1, 1) #shape:  (1348, 1, 1) ; 1348 tetrehedral element number in cylinder.mesh file
    print "shape: ",shape
    out = {
        'fmax' : vf * nm.tile(fmax, shape), # 1348 elemanı ve değeri 10.0 olan bir vector yaratır.
        'eps_opt' : nm.tile(eps_opt, shape), # 1348 elemanı ve değeri 0.01 olan bir vector yaratır.
        's' : nm.tile(s, shape),  # 1348 elemanı ve değeri 1.0 olan bir vector yaratır.
        'fdir' : nm.tile(fdir, shape), # 1348 satırlı değeri [1,0,0] olan bir array yaratır.
        'act' : nm.tile(act, shape), # 1348 elemanı ve değeri act olan bir vector yaratır.
    }

    return out
"""
#######################################################################################################
"""
.  mode is the actual evaluation mode, default is ‘eval’;

.  **kwargs: keyword arguments.

.  which ve vf arguments'ları extra argumentlerdir. Bu yüzden 'get_pars_fibres' fonksiyonunu çağırırken lambda kullanılır.

.  def get_pars_fibres fonksiyonun döndürdüğü değer get_pars_fibres1 ve get_pars_fibres2 olmak üzere iki değere atandı.
"""
"""
    'get_pars_fibres1' : (lambda ts, coors, mode=None, **kwargs:
                          get_pars_fibres(ts, coors, mode=mode, which=0,
                                          vf=vf_fibres1, **kwargs),),
"""

functions = {

    'get_bc' : (lambda ts, coors, mode=None, **kwargs:
                get_boundary(ts,coors, mode=mode),),
}
#######################################################################################################
"""
.  u : -> Unknown field, name of which is displacement, 
       -> 0; means order in the global vector of unknowns,. That is, it’s associated DOF name is 'u.0'
     
.  v : -> test variable of the weak formulation is called 'v'.
       -> a test variable must specify the unknown(u) it corresponds to.
"""

variables = {
    'u' : ('unknown field', 'displacement', 0),
    'v' : ('test field', 'displacement', 'u'),
}

#######################################################################################################
"""
.  Omega : all the elements in the mesh. In other words, the elemental domain over which the PDE is solved 

.  Left and Right : -> define surfaces upon which the boundary conditions will be applied.
                    -> Because of D = 3, facets are faces.(D-1) 
                    -> vertices in (expr) means faces are specified wrt x values and they are left and right faces of cylinder.

"""
regions = {
    'Omega' : 'all',
    'Left' : ('vertices in (x < -5.9)', 'facet'),
    'Right' : ('vertices in (x > 5.9)', 'facet'),
    'point15': ('vertex 14', 'vertex'),
}

#######################################################################################################

# Dirichlet BC.
"""
.  l : the displacements(u and v) of the nodes in the Left region are prevented or set to zero 

.  Notice that; if we write {'u.0' : 0.0}, horizontal(x) displacements of the nodes are prevented.
                if we write {'u.1' : 0.0}, vertical(y) displacements of the nodes are prevented.
                if we write {'u.0': 'get_bc' , 'u.[1,2]' : 0.0}, predisplacement for x displacement from get_bc function is applied,
                                                                 y and z displacements are prevented.

"""

ebcs = {
    'l' : ('Left', {'u.0' : 0.0}),
    'r' : ('Right', {'u.0': 'get_bc' }),
    'p_15' : ('point15', {'u.all' : 0.0}),
}


    

#######################################################################################################

# Balance of forces.
"""
.  Here we are using a 1nd order quadrature over a 3 dimensional space.

"""
integral_1 = {
    'name' : 'i',
    'order' : 1,
}

#######################################################################################################
"""
           + dw_tl_fib_a.i.Omega( f1.fmax, f1.eps_opt, f1.s, f1.fdir, f1.act,
                                  v, u )
"""
equations = {
    'balance'
        : """dw_tl_he_neohook.i.Omega( solid.mu, v, u )
           + dw_tl_bulk_penalty.i.Omega( solid.K, v, u )
           = 0""",
           

}

#######################################################################################################
"""   
    This will be called after the problem is solved.

    Parameters
    ----------
    out : dict
        The output dictionary, where this function will store additional
        data.
    problem : Problem instance
        The current Problem instance.
    state : State instance
        The computed state, containing FE coefficients of all the unknown
        variables.
    extend : bool
        The flag indicating whether to extend the output data to the whole
        domain. It can be ignored if the problem is solved on the whole
        domain already.

    Returns
    -------
    out : dict
        The updated output dictionary.
    

.  ‘el_avg’ : -> The element average mode results in an array of a quantity averaged in each element of a region. 
              -> This is the mode for postprocessing.

"""
"""
def Factivation(t):
  t_d = 0.0
  t_a = 0.4
  q = 3.5
  t_p = 0.8
  S = 200
  fa = (1+nm.sin(nm.pi*(t-t_d)/t_a - nm.pi/2))/2
  if t <= t_d:
     f_act = 0.0
    
  if t_d < t <= t_d + t_a:
     f_act = pow(fa,q) 
       
  if t_d + t_a < t <= t_p:
     f_act = 1.0
       
  if t_p < t:
     f_act = nm.exp(-S*(t-t_d-t_p))
  return f_act
"""

# Boundary Condition value kademeli olarak uygulanır. Mesela 5 ise 1,2,3,4,5 şeklinde 5 step'te uygulanır.

def get_boundary(ts,coors, mode=None):
    if bc_value >= 0:
        if ts.time*160 <= (bc_value):
            a = ts.time * 160
        else:
            a = bc_value
            
    else:
        if ts.time*(-160) >= (bc_value):
            a = ts.time * (-160)
        else:
            a = bc_value
    
    return a


def stress_strain(out, problem, state, extend=False):
    from sfepy.base.base import Struct

    ev = problem.evaluate 
                          
    
    strain1 = ev('dw_tl_he_neohook.i.Omega( solid.mu, v, u )',
                mode='el_avg', term_mode='strain')
    out['green_strain'] = Struct(name='output_data',
                                 mode='cell', data=strain1, dofs=None)

    stress1 = ev('dw_tl_he_neohook.i.Omega( solid.mu, v, u )',
                mode='el_avg', term_mode='stress')
    out['neohook_stress'] = Struct(name='output_data',
                                   mode='cell', data=stress1, dofs=None )

    stress2 = ev('dw_tl_bulk_penalty.i.Omega( solid.K, v, u )',
                mode='el_avg', term_mode= 'stress')
    out['bulk_stress'] = Struct(name='output_data',
                                mode='cell', data=stress2, dofs=None)
    """                            
    stress3 = ev('dw_tl_fib_a.i.Omega( f1.fmax, f1.eps_opt, f1.s, f1.fdir, f1.act, v, u)',
                mode='el_avg', term_mode= 'stress')
    out['active_stress'] = Struct(name='output_data',
                                mode='cell', data=stress3, dofs=None)
  
    out['Piola_Kirchhoff_stress'] = Struct(name='output_data',
                                mode='cell', data=stress1 + stress3, dofs=None)
                                
    out['Total_stress'] = Struct(name='output_data',
                                mode='cell', data=stress1 + stress2 + stress3, dofs=None)
    """
    return out

#######################################################################################################
# Solvers etc.
"""
name : str

   The name referred to in problem description options.

kind : str

   The solver kind, as given by the name class attribute of the Solver subclasses.

"""

solver_0 = {
    'name' : 'ls',
    'kind' : 'ls.scipy_direct', # Direct sparse solver from SciPy.
}

#######################################################################################################
"""
class sfepy.solvers.nls.Newton

.  Solves a nonlinear system f(x) = 0 using the Newton method with backtracking line-search, 
starting with an initial guess x^0.

.  Kind: ‘nls.newton’

.  Parameters
     'i_max'      : The maximum number of iterations.
     'eps_a'      : The absolute tolerance for the residual, i.e. ||f(x^i)||.
     'eps_r'      : The relative tolerance for the residual, i.e. ||f(x^i)|| / ||f(x^0)||.
     'macheps'    : The float considered to be machine “zero”.
     'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
     'ls_red'     : The linear system solution error should be smaller than (eps_a * lin_red) 
     'ls_red_warp': The step reduction factor in case of failed residual assembling.
     'ls_on'      : Start the backtracking line-search by reducing the step, if ||f(x^i)|| / ||f(x^{i-1})|| is larger than ls_on.
     'ls_min'     : The minimum step reduction factor
     'check'      : If >= 1, check the tangent matrix using finite differences. If 2, plot the resulting sparsity patterns.
     'delta'      : If check >= 1, the finite difference matrix is taken as Aij = ....
"""


solver_1 = {
    'name' : 'newton',
    'kind' : 'nls.newton',

    'i_max'      : 7,
    'eps_a'      : 1e-7,
    'eps_r'      : 1.0,
    'macheps'    : 1e-16,
    'lin_red'    : 1e-2, # Linear system error < (eps_a * lin_red).
    'ls_red'     : 0.1,
    'ls_red_warp': 0.001,
    'ls_on'      : 1.1,
    'ls_min'     : 1e-5,
    'check'      : 0,
    'delta'      : 1e-6,
}
#######################################################################################################
"""
class sfepy.solvers.ts_solvers.SimpleTimeSteppingSolver
 
.  Implicit time stepping solver with a fixed time step.

.  Kind: ‘ts.simple’

.  Parameters:       
     n_step : Iteration number -> iteration number starts 0, so iteration stops n_step = 20
     t0     : The initial time.
     t1     : The final time.
     dt     : Time-Step       -> Used if n_step is not given. When dt=None indicated like below, dt = (t1-t0)/(n_step - 1) is calculated in ts.py file.
     time   : CurrentTime   -> time += dt, Iteration is over when time = t1
     nt     : Normalized Time  -> nt += (time-t0)/(t1 - t0)
"""
solver_2 = {
    'name' : 'ts',
    'kind' : 'ts.simple',

    't0'    : 0,
    't1'    : 1,
    'dt'    : None,
    'n_step' : 161, # has precedence over dt!
}
