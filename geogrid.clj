(ns
    geogrid
  "A geographic grid protocol and associated function"
  (:use geoprim)
  (:require [thi.ng.geom
             [core                :as geom]
             [matrix              :as matrix]]
            [thi.ng.geom.viz.core :as viz]
            [thi.ng.geom.svg.core :as svg]
            [thi.ng.math.core     :as math]
            [thi.ng.ndarray.core  :as ndarray]
            [thi.ng.color.core    :as col]))

;;(set! *warn-on-reflection* true)
;;(set! *unchecked-math* :warn-on-boxed)

;; GRID
(defprotocol
    grid
  "A grid of geographic data - evenly spaced in lat/lon"
  (dimension-pix
    [x]
    "Return the [width height] of a region in pixels")
  (eassou-res
    [x]
    "Returns [eas-res sou-res] of the grid")
  (corner
    [x]
    "Returns the grid's start position (norwes corner) as a `geopoint`")
  (data
    [x]
    "Return a vector of the grid data")
  (subregion
    [x
     region]
    "Return a subregion of the grid - consisting of all points in the region bounds"))

(defn
  params
  "Extract all the parameters into a list
  [width-pix
   height-pix
   eas-res
   sou-res
   norwes-point]
  Used to initialized a new grid in for example `geogrid4seq/build-grid`"
  [given-grid]
  (let [[width
         height]     (-> given-grid
                         dimension-pix)
        [eas-res
         sou-res]    (-> given-grid
                         dimension-pix)
        norwes-point (-> given-grid
                         corner)]
    [width
     height
     eas-res
     sou-res
     norwes-point]))

(defn
  point-to-pix
  "Given `grid` 
  Return the pixel within which a given `point` lies.
  .[Note]
  * pixel coordinates are zero indexed
  * coordinates can be negative
  * the first negative pixel position is `-1`
  * the value is a fractional `double`"
  [given-point
   given-grid]
  (let[[^double
        eas-point
        ^double
        sou-point]  (as-eassou
                      given-point)
       [^double
        eas-corner
        ^double
        sou-corner] (as-eassou
                      (->
                        given-grid
                        corner))
       [^double
        eas-res
        ^double
        sou-res]    (->
                      given-grid
                      eassou-res)]
    [(->
       eas-point
       (-
         eas-corner)
       (/
         eas-res))
     (->
       sou-point
       (-
         sou-corner)
       (/
         sou-res))]))

(defn
  adjusted-crop-region-to-grid
  "Adjust a crop region so that it's aligned with the image and
  such that the sides are an even multiple of the pixels.
  Adjustment will always makes the crop region a bit larger
  At most by  a fraction of a pixel (unless it already fits)
  And it will always include the full requested area.
  == Return:
  `:crop-region-pixel-offsets`::
  with values `:norwes` and `:soueas` which are pixel offsets
  of the north west and south east corners of the region after adjustment
  `:overruns`::
  how far the new crop region overruns the requested area
  This is units of fractions-of-a-pixel (ie. it's always less than `1.0`)"
  [crop-region
   given-grid]
  (let[[norwes
        noreas
        soueas
        souwes]       (four-corners
                        crop-region)
       norwes-pix     (->
                        norwes
                        (point-to-pix
                          given-grid))
       norwes-rounded (->>
                        norwes-pix
                        (mapv
                          #(->
                             ^double %
                             Math/floor
                             int)))
       norwes-overrun (mapv
                        -
                        norwes-pix
                        norwes-rounded)
       soueas-pix     (point-to-pix
                        soueas
                        given-grid)
       soueas-rounded (->>
                        soueas-pix
                        (mapv
                          #(->
                             ^double
                             %
                             dec
                             Math/ceil
                             long)))
       soueas-overrun (mapv
                        -
                        (mapv
                          inc
                          soueas-rounded)
                        soueas-pix)]
    {:crop-region-pixel-offsets {:start-x (first
                                            norwes-rounded)
                                 :ended-x (first
                                            soueas-rounded)
                                 :start-y (second
                                            norwes-rounded)
                                 :ended-y (second
                                            soueas-rounded)}
     :overruns {:top    (->
                          norwes-overrun
                          second)
                :bottom (->
                          soueas-overrun
                          second)
                :right  (->
                          soueas-overrun
                          first)
                :left   (->
                          norwes-overrun
                          first)}}))
;;.[EXAMPLE]
#_ 
(let[minnan-region (region
                     (point
                       26.23
                       116.47)
                     (point
                       21.7
                       125))
     rain-grid     (->
                     simple-grid ;; TODO FIX
                     100
                     100
                     1.0
                     1.0
                     (point
                       25
                       120)
                     [])]
  (adjusted-crop-region-in-image-coords
    minnan-region
    rain-grid))
;;
;; => {:crop-region-pixel-offsets {:norwes [-4 -2], :soueas [4 3]},
;;     :overruns
;;     {:top 0.769999999999996,
;;      :bottom 0.7000000000000028,
;;      :right 0.0,
;;      :left 0.4700000000000273}}

(defn
  normalized-data
  "calls `data` and normalizes the values
  .. to the `0.0-1.0` range"
  ([grid
    magnitude]
   (mapv
     #(/
        %
        magnitude)
     (data
       grid)))
  ([grid]
   (let[data-seq  (data
                    grid)
        max-value (apply
                    max
                    data-seq)
        min-value (apply
                    min
                    data-seq)]
     (normalized-data
       grid
       (max
         (Math/abs
           max-value)
         (Math/abs
           min-value))))))

(defn
  covered-region
  "Uses the region parameters to calculate
  and then return the region covered by the grid"
  [input-grid]
  (let [norwes-corner (-> input-grid
                          corner)
        [width
         height]      (-> input-grid
                          dimension-pix)
        [eas-res
         sou-res]     (-> input-grid
                          eassou-res)]
      (let[soueas-corner (point-eassou (+ (-> norwes-corner
                                              :eas)
                                          (* width
                                             eas-res))
                                       (+ (-> norwes-corner
                                              :sou)
                                          (* height
                                             sou-res)))]
    (region norwes-corner
            soueas-corner))))
