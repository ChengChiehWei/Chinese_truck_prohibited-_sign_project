
Open "main.m" with matlab,then hit "run" button.

This function will find targets in pictures under folder "database" by comparing templates in folder "template".
A matched target will be marked with blue square.
However, with very limited time, I adjusted parameters for "template1.jpg" (the standard Chinese "truck prohibited" sign) only.
Might not work very well with other templates.

Since my algorithm adopts optimization methods, you might have to try couple of times to see it locks on right targets.
(Because the optimization methods might fall into local minimum sometime. But generally speaking, it still has a higher probability to lock on right targets.)

