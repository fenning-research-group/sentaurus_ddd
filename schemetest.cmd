; file used to test scheme language (run through Sentaurus structure editor)
; command: sde -e -l schemetest.cmd

(define dopingmodel "uniform")
(define dopingmodel "erf")

(cond 
  ((string=? dopingmodel "uniform")
	(begin
		(newline)
		(display "***********************************************************")
		(newline)
		(display "Using an uniform doping profile")
		(newline)
		(display "***********************************************************")
		(newline)
	)
  )
  ((string=? dopingmodel "erf")
	(begin
		(display "Using an error function doping profile")
	)
  )
)