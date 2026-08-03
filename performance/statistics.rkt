#lang racket/base

(require racket/match)

(provide (struct-out stats)
         make-stats
         statistical-analysis)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Simple statistical analysis tools

;; Represents random variable with standard deviation and basic
;; information about source data set
(struct stats
  (avg    ;; Arithmetic average
   stdev  ;; Standard deviation
   min    ;; Minimum value
   max    ;; Maximum value
   count  ;; Number of samples
   )
  #:methods gen:custom-write
  ((define (write-proc v p m)
     (match-define (stats avg stdev min max count) v)
     (display
      (format "~a ±~a [~a] (~a .. ~a)"
              (real->decimal-string avg 6)
              (real->decimal-string stdev 6)
              count
              min
              max)
      p))))

;; Converts a list of numbers into the statistical information about
;; that list
(define (make-stats lst)
  (define cnt (length lst))
  (define sum (foldl + 0 lst))
  (define avg (/ sum cnt))
  (define minv (apply min lst))
  (define maxv (apply max lst))
  (define variance
    (/ (for/sum ((v (in-list lst)))
         (expt (- v avg) 2))
       cnt))
  (define stdev (sqrt variance))
  (stats avg stdev minv maxv cnt))

;; Returns the list with one value outside of given sigma-multiple
;; range from the list. The farthest value is removed. Returns a new
;; list and an alement that was removed or #f.
(define (remove-the-most-outside-sigma raw-values #:max-sigma (max-sigma 3))
  (define raw-stats (make-stats raw-values))
  (define max-delta (* max-sigma (stats-stdev raw-stats)))
  (define raw-average (stats-avg raw-stats))
  (for/fold ((clean-values '())
             (current-removed #f)
             (current-delta 0)
             #:result (values clean-values
                              current-removed))
            ((val (in-list raw-values)))
    (define val-delta (abs (- val raw-average)))
    (if (< val-delta max-delta)
        (values (cons val clean-values)
                current-removed
                current-delta)
        (if (and current-removed
                 (< current-delta val-delta))
            (values (cons current-removed clean-values)
                    val
                    val-delta)
            (values (cons val clean-values)
                    current-removed
                    current-delta)))))

;; Performs statistical analysis of raw input values and returns a new
;; list of input values with extremes removed and the statistical
;; analysis of remaining values.
(define (statistical-analysis raw-values #:max-sigma (max-sigma 3))
  (define clean-values
    (let loop ((last-values raw-values))
      (define-values (new-values removed)
        (remove-the-most-outside-sigma last-values #:max-sigma max-sigma))
      (cond (removed
             (loop new-values))
            (else
             new-values))))
  (define clean-stats (make-stats clean-values))
  (values clean-values
          clean-stats))
