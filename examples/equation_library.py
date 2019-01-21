# Define the equations of motion
euler_1d ="""
# Primary Equations of motion here
ddt(:rho:)  =  -ddx(:rho:*:u:)
ddt(:rhou:) =  -ddx(:rhou:*:u: + :p: - :tau:)
ddt(:Et:)   =  -ddx( (:Et: + :p: - :tau:)*:u: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u:) ) * ( :gamma: - 1.0 )
# Artificial bulk viscosity (old school way)
:div:       =  ddx(:u:) 
:beta:      =  gbar(abs(ring(:div:))) * :rho: * 7.0e-2
:tau:       =  :beta:*:div:
"""


euler_3d ="""
# Primary Equations of motion here
ddt(:rho:)  =  -div(:rho:*:u:  ,  :rho:*:v: , :rho:*:w:)
ddt(:rhou:) =  -div(:rhou:*:u: - :tauxx: , :rhou:*:v: - :tauxy:, :rhou:*:w: - :tauxz:)
ddt(:rhov:) =  -div(:rhov:*:u: - :tauxy: , :rhov:*:v: - :tauyy:, :rhov:*:w: - :tauyz:)
ddt(:rhow:) =  -div(:rhow:*:u: - :tauxz: , :rhow:*:v: - :tauyz:, :rhow:*:w: - :tauzz:)
ddt(:Et:)   =  -div( (:Et: - :tauxx:)*:u: - :tauxy:*:v: - :tauxz:*:w: , (:Et: - :tauyy:)*:v: -:tauxy:*:u: - :tauyz:*:w: ,  (:Et: - :tauzz:)*:w: - :tauxz:*:u: - :tauyz:*:v: )
# Conservative filter of the EoM
:rho:       =  fbar( :rho:  )
:rhou:      =  fbar( :rhou: )
:rhov:      =  fbar( :rhov: )
:rhow:      =  fbar( :rhow: )
:Et:        =  fbar( :Et:   )
# Update the primatives and enforce the EOS
:u:         =  :rhou: / :rho:
:v:         =  :rhov: / :rho:
:w:         =  :rhow: / :rho:
:p:         =  ( :Et: - .5*:rho:*(:u:*:u: + :v:*:v: + :w:*:w:) ) * ( :gamma: - 1.0 )
# Divergence and cross derivatives
:div:       =  div(:u:,:v:,:w:) 
[:ux:,:uy:,:uz:] = grad(:u:)
[:vx:,:vy:,:vz:] = grad(:v:)
[:wx:,:wy:,:wz:] = grad(:w:)
# Artificial bulk viscosity 
:mu:        =  gbar( ringV(:u:,:v:,:w:) * :rho: ) * 1.0e-1 + mu0
:beta:      =  gbar( abs(ring(:div:)) * :rho: )  * 7.0e-3
:taudia:    =  (:beta:-2./3.*:mu:) *:div: - :p:
:tauxx:     =  2.0*:mu:*:ux:   + :taudia:
:tauyy:     =  2.0*:mu:*:vy:   + :taudia:
:tauzz:     =  2.0*:mu:*:wz:   + :taudia:
:tauxy:     = :mu:*(:uy:+:vx:)
:tauxz:     = :mu:*(:uz:+:wx:)
:tauyz:     = :mu:*(:vz:+:wz:)
:cs:  = sqrt( :p: / :rho: * :gamma: )
:dt: = dt.courant(:u:,:v:,:w:,:cs:)*1.0
:dt: = numpy.minimum(:dt:,0.2 * dt.diff(:beta:,:rho:))
:dt: = numpy.minimum(:dt:,0.2 * dt.diff(:mu:,:rho:))
"""
