// Check if we've already defined a MIP consistency header file (don't redefine)
#ifndef MIPCONSISTENCY_H
#define MIPCONSISTENCY_H

// Include statements go here


// Now declare the function
bool IsMIP(recob::track &track, art::Event & evt);

#endif
