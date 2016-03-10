#include "sampling_event.h"
#include "Sample_With_Replacement.h"
#include "Sample_without_Replacement_1Pass.h"
// #include "Sample_without_Replacement.h" // For now, this isn't quite ready yet

// TODO: right now have no way of returning NULL if error?

SamplingEvent *SamplingEvent::create_SamplingEvent(SamplingInput *sampling_input, curandGenerator_t gen)
	{
	if (sampling_input->sampling_scheme == 1)
		return new Sample_With_Replacement(sampling_input, gen);

	if (sampling_input->sampling_scheme == 2)
		return new Sample_without_Replacement_1Pass(sampling_input, gen);

 /*
// not ready 
  if (sampling_input->sampling_scheme == 3)
    return new Sample_without_Replacement(sampling_input, gen);

*/
	}


