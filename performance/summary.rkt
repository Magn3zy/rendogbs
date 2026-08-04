#lang racket/base

(require racket/cmdline
         racket/path
         racket/function
         racket/string
         racket/port
         racket/match
         "statistics.rkt"
         racket/list
         plot/pict
         pict
         racket/format)

(current-directory "ddr4")

;(define dirnames (command-line #:args paths paths))

(define names '(pp csa td pm hv sc ta))
(define dirnames '("results_pp/" "results_csa/" "results_td/" "results_pm/" "results_hv/" "results_sc/" "results_ta/"))
(define sizes '(230440408 3098728024 10636951312 845112719 4278575539 6786457495 14762960853))

(define (parse-result fname)
  (for/hash ((a (in-list
                 (filter (λ (v)
                           (and (not (null? v))
                                (not (null? (cdr v)))))
                         (map (curryr string-split ":")
                              (with-input-from-file fname port->lines))))))
    (values (car a) (string-join (cdr a) ":"))))

(define (parse-data fname)
  (define lst (map string-split (with-input-from-file fname port->lines)))
  (define head
    (let ((h0 (car lst)))
      (cons (substring (car h0) 1)
            (cdr h0))))
  (define body
    (for/fold ((lt 0)
               (res #f)
               #:result res)
              ((row (in-list (cdr lst)))
               #:when (and (= (length row)
                              (length head))
                           (>= (string->number (car (reverse row))) lt)))
      (values (string->number (car (reverse row)))
              row)))
  #;(define body
    (for/last ((b0 (in-list (cdr lst)))
               #:when (= (length b0) (length head)))
      b0))
  (define baseline
    (for/first ((b0 (in-list (cdr lst)))
               #:when (= (length b0) (length head)))
      b0))
  (values (for/hash ((k (in-list head))
                     (v (in-list body)))
            (values k (string->number v)))
          (for/hash ((k (in-list head))
                     (v (in-list baseline)))
            (values k (string->number v)))))

(define (process-benchmark-results dirname)
  (define runs
    (for/list ((rname (in-directory dirname (lambda (_) #f)))
               #:when (and (path-has-extension? rname ".result")
                           (file-exists? (path-replace-extension rname #".data"))))
      (define dname (path-replace-extension rname #".data"))
      (define-values (_1 bname _2) (split-path (path-replace-extension rname #"")))
      (define result (parse-result rname))
      (define-values (data baseline) (parse-data dname))
      (list bname result data baseline)))
  (define sname (build-path dirname "runs.txt"))
  (define onlyruns
    (with-output-to-file sname #:exists 'replace
      (thunk
       (displayln "#DURATION THREADS MEMORY USERMS SYSTEMMS TOTALMS")
       (for/list((run (in-list runs)))
         (match-define (list bname result data baseline) run)
         (define duration (- (hash-ref data "TS") (hash-ref baseline "TS") 2)) ; 2s from --benchmark-sleep
         (define num-threads (- (hash-ref data "MAXTHREADS") (hash-ref baseline "MAXTHREADS")))
         (define memory (- (hash-ref data "MAXMEMORY") (hash-ref baseline "MAXMEMORY")))
         (define user/ms (- (hash-ref data "USER") (hash-ref baseline "USER")))
         (define system/ms (- (hash-ref data "SYSTEM") (hash-ref baseline "SYSTEM")))
         (define total/ms (- (hash-ref data "TOTAL") (hash-ref baseline "TOTAL")))
         ;(displayln (format "# ~a: ~a" bname (hash-ref result "Ran with options" "COMMAND-LINE ARGUMENTS MISSING")))
         (displayln (format "~a ~a ~a ~a ~a ~a" duration num-threads memory user/ms system/ms total/ms))
         (list duration num-threads memory user/ms system/ms total/ms)))))
  (list (map first onlyruns) ; duration
        (map second onlyruns) ; threads
        (map third onlyruns) ; memory
        (map sixth onlyruns) ; total[ms]
        ))

(define results
  (for/list ((dirname (in-list dirnames)))
    (process-benchmark-results dirname)))

(define-values (tstats mstats pstats cstats)
  (for/lists (a b c d)
             ((res (in-list results)))
    (values (make-stats (first res))
            (make-stats (map (λ (v) (/ v 1024 1024 1024)) (third res)))
            (make-stats (map (λ (v) (/ v 1)) (second res)))
            (make-stats (map (λ (v) (/ v 1000)) (fourth res))))))

(hc-append
 (plot
  `(
    ,(tick-grid)
    ,(points
      (for/list ((rts (in-list tstats))
                 (x (in-list sizes)))
        (list x (stats-avg rts)))
      #:color "red"
      #:fill-color "magenta"
      #:sym 'fullcircle)
    ,@(for/list ((label (in-list names))
                 (x (in-list sizes))
                 (rts (in-list tstats)))
        (point-label (vector x (stats-avg rts)) (format "~a" label)))
    ,(error-bars
      (for/list ((rts (in-list tstats))
                 (x (in-list sizes)))
        (list x (stats-avg rts) (* 2 (stats-stdev rts))))))
  #:x-min 0
  #:y-max (* 1.05 (apply max (map stats-max tstats)))
  #:x-max (* 1.1 (apply max sizes))
  #:x-label "Genome Size"
  #:y-label "Running Time [s]"
  #:y-min 0)

 (plot
  `(
    ,(tick-grid)
    ,(points
      (for/list ((rts (in-list mstats))
                 (x (in-list sizes)))
        (list x (stats-avg rts)))
      #:color "red"
      #:fill-color "magenta"
      #:sym 'fullcircle)
    ,@(for/list ((label (in-list names))
                 (x (in-list sizes))
                 (rts (in-list mstats)))
        (point-label (vector x (stats-avg rts)) (format "~a" label)))
    ,(error-bars
      (for/list ((rts (in-list mstats))
                 (x (in-list sizes)))
        (list x (stats-avg rts) (* 2 (stats-stdev rts))))))
  #:x-min 0
  #:y-max (* 1.05 (apply max (map stats-avg mstats)))
  #:x-max (* 1.1 (apply max sizes))
  #:x-label "Genome Size"
  #:y-label "Memory Usage [GB]"
  #:y-min 0))

(hc-append
 (plot
  `(
    ,(tick-grid)
    ,(points
      (for/list ((rts (in-list cstats))
                 (x (in-list sizes)))
        (list x (stats-avg rts)))
      #:color "red"
      #:fill-color "magenta"
      #:sym 'fullcircle)
    ,@(for/list ((label (in-list names))
                 (x (in-list sizes))
                 (rts (in-list cstats)))
        (point-label (vector x (stats-avg rts)) (format "~a" label)))
    ,(error-bars
      (for/list ((rts (in-list cstats))
                 (x (in-list sizes)))
        (list x (stats-avg rts) (* 2 (stats-stdev rts))))))
  #:x-min 0
  #:y-max (* 1.05 (apply max (map stats-avg cstats)))
  #:x-max (* 1.1 (apply max sizes))
  #:x-label "Genome Size"
  #:y-label "CPU Time [s]"
  #:y-min 0)

 (plot
  (list
   (tick-grid)
   (points
    (for/list ((rts (in-list pstats))
               (x (in-list sizes)))
      (list x (stats-max rts)))
    #:color "red"
    #:fill-color "magenta"
    #:sym 'fullcircle)
   #;(error-bars
      (for/list ((rts (in-list pstats))
                 (x (in-list sizes)))
        (list x (stats-avg rts) (* 2 (stats-stdev rts))))))
  #:x-min 0
  #:y-max (* 1.11 (apply max (map stats-avg pstats)))
  #:x-max (* 1.05 (apply max sizes))
  #:x-label "Genome Size"
  #:y-label "Threads"
  #:y-min 0))

(displayln "#SIZE THREADS TIME TIMESTDEV MEMORY MEMSTDEV CPUTM CPUTMSTDEV LABEL")

(define rows
  (sort
   (for/list ((ts (in-list tstats))
              (ms (in-list mstats))
              (ps (in-list pstats))
              (cs (in-list cstats))
              (label (in-list names))
              (x (in-list sizes)))
     (list
      x
      (stats-max ps)
      (real->decimal-string (stats-avg ts)) (stats-stdev ts)
      (real->decimal-string (stats-avg ms)) (stats-stdev ms)
      (real->decimal-string (stats-avg cs)) (stats-stdev cs)
      label
      ))
   (lambda (a b)
     (< (car a) (car b)))))

(for ((row (in-list rows)))
  (displayln (string-join (map ~a row) " ")))
