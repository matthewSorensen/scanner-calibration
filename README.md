# Scanner Calibration


Although  laser scanner systems' f-theta lenses are designed to have a linear relationship between beam
angle and spot location, essentially all high precision applications require additional empirical calibration
to minimize optical distortion across the full field of view. This repository contains scripts to generate
test patterns, analyze high-resolution scans of the test patterns, and generate distortion compensation models
for on-the-fly error correction of scanner trajectories. 

## Warning

Powerful lasers can, and will, cause instantaneous damage to skin, eyes, and other organic tissues. This document
is not intended to instruct users on laser safety, and users must educate themselves on the risks and mitigations
before operating any laser equipment.

## Supported Scanner Control Systems

At the moment, this repository only supports [pewpew-laser](https://github.com/matthewSorensen/pewpew-laser), an
unreleased scanner and laser control system running on the teensy 4.x microcontroller series.

## Supplies and Equipment

The following supplies and equipment are required to complete the lens calibration process:

  * A galvo-scanned laser system running a compatible control system, and suitable safety equipment for its
      operation.
  * A dimensionally stable material capable of high-resolution, high-contrast marking. Low-cost
      black cardstock manufactured with a white core is a reasonable material for 1064nm laser calibration.
      Black anodized aluminum is an ideal material, but is significantly more expensive, particularly for larger
      sheets. Cyanotype paper may be suitable for UV and violet laser sources.
  * A flatbed scanner capable of imaging the laser-marked test pattern with minimal distortion.

