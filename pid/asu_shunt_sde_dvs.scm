
(sde:clear)

(sde:set-process-up-direction "+z")

(define MNAME       "Metal"             )
(define MBULK       "Silicon"	        )
(define MCAP        "Nitride"	        )
(define PARAM       "MetalConductivity" )
(define DEPTH       20.0                )
(define INIT	    "a=1.2"             )
(define FUNCTION    "1.0-a^(x-20.0)"    )
(define DEFAULT     0.0                 )

;"a=1.2" "1.0-a^(x-20.0)" 0.0


(sdegeo:create-rectangle
	(position 0 0 0)
	(position DEPTH 4 0)
	MBULK
	"Bulk")

(sdegeo:create-rectangle
	(position -2 0 0)
	(position 0 4 0)
	MCAP
	"Cap")

(sdegeo:create-rectangle
	(position  0 1.5 0)
	(position DEPTH 2.5 0)
	MNAME
	"MetalPlug")

#if "@mode@" == "function"
(sdedr:define-refeval-window
	"Window.Sigma"
	"Rectangle"
	(position 0 1.5 0)
	(position DEPTH 2.5 0) )

(sdedr:define-analytical-profile
	"Prof.Sigma"
	"MetalConductivity"
	INIT FUNCTION DEFAULT "general")

(sdedr:define-analytical-profile-placement
	"Place.Sigma"
	"Prof.Sigma"
	"Window.Sigma"
	"Positive"
	"Replace"
	"Eval"
	"Window.Sigma"
	0.0
	"evalwin")
#endif

#if "@mode@" == "1d_profile"
(sdedr:define-1d-external-profile
	"Prof.Sigma"
	"conductivity.plx"
	"Range" 0 DEPTH
	"Erf" "Length" 0)

(sdedr:define-refeval-window
	"BaseLine.Sigma"
	"Line"
	(position 0.0 1.5 0)
	(position 0.0 2.5 0) )

(sdedr:define-analytical-profile-placement
	"Place.Sigma"
	"Prof.Sigma"
	"BaseLine.Sigma"
	"Negative"
	"Replace"
	"Eval"
	""
	0.0
	"evalwin")
#endif

(sdedr:define-refinement-size
	"RS.Global"
	0.1 0.1 0.0
	0.01 0.01 0.0)

(sdedr:define-refinement-placement
	"RP.Global"
	"RS.Global"
	"RefEvalWin.Global")

(sdedr:define-refinement-size
	(string-append "RS." MNAME)
	0.1  0.1  0.0
	0.01 0.01 0.0)

(sdedr:define-refinement-material
	(string-append "RP." MNAME)
	(string-append "RS." MNAME)
	MNAME)

(sdedr:define-refinement-material
	(string-append "RP." "Silicon")
	(string-append "RS." MNAME)
	"Silicon")

(sde:build-mesh "" "n@node@")
