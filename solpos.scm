#!r6rs
;; A solar position calculator based on the 6S library
;; (See http://6s.ltdri.org/index.html).

;; The library exports just one function: solar-position, which given
;; a SRFI 19 date and a ground position in lat/long, returns the
;; zenithal and azimuthal angles for the solar position in sexagesimal
;; degrees.

(library (solpos)
         (export solar-position)
         (import (rnrs)
                 (srfi :19)
                 (rnrs base (6)))

(define pi (* (atan 1.0) 4.0))

(define-syntax trig-coeffs
  (lambda (stx)
    (define (alternate stx-angle stx-coeff)
      (datum->syntax stx-angle
                     (do ((i 2 (+ i 1))
                          (coef (syntax->datum stx-coeff) (cdr coef))
                          (result '() (cons
                                       (list '* (car coef)
                                             (list (if (odd? i) 'sin 'cos)
                                                   (list '*
                                                         (floor (/ i 2))
                                                         (syntax->datum stx-angle))))
                                            result)))
                         ((null? coef) (reverse result)))))
    (syntax-case stx (:)
      [(_ angle : zero . coeff)
       #`(+ zero #,@(alternate #'angle #'coeff))])))

(define (solar-position date lon lat)
  (define doy (+ 0 (date-year-day date)))
  (define decimal-hour (+ (date-hour date)
                          (/ (date-minute date)
                             60)))
  (define fac (/ pi 180))
  (define tsm (+ decimal-hour
                 (/ lon 15)))
  (define xla (* fac lat))
  (define tet (/ (* pi 2 doy )
                 365.0))
  (define sin-tet (sin tet))
  (define cos-tet (cos tet))
  
  (define et (/ (* (trig-coeffs tet : 0.0000750 0.0018680 -0.0320770
                                -0.0146150 -0.0408490)
                   12 60)
                pi))
  
  (define tsv (+ tsm (/ et 60) -12))
  (define ah (* tsv 15 fac))
  
  (define delta (trig-coeffs tet : 0.0069180 -0.3999120 0.0702570
                             -0.0067580 0.0009070 -0.0026970 0.0014800))
  
  (define amuzero (+ (* (sin xla)
                        (sin delta))
                     (* (cos xla)
                        (cos delta)
                        (cos ah))))
  
  (define elev (asin amuzero))
  (define pre-az (/ (* (cos delta)
                       (sin ah))
                    (cos elev)))
  (define az (cond [(<= pre-az -1.0) -1.0]
                   [(> pre-az 1.0) 1.0]
                   [#t pre-az]))
  
  
  (define caz (+ (* (cos xla)
                    (sin delta)
                    -1)
                 (/ (* (sin xla)
                       (cos delta))
                    (cos elev))))

  (define azim (asin az))
  (define two-times-pi (* pi 2))

  (if (<= caz 0)
      (set! azim (- pi azim))
                                        ;else
      (when (<= az 0)
        (set! azim (+ azim two-times-pi))))

  (set! azim (+ azim pi))
  (when (> azim two-times-pi)
    (set! azim (- azim two-times-pi)))


  (values (- 90.0 (/ elev fac))
          (/ azim fac)))

)
