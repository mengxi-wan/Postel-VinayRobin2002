#----------------------------------------------------------------
# Optimzation requires the parameter argument to be vector
# This function converts the input parameter vector into a struct 
# to generate worker panel and simulated moments
#----------------------------------------------------------------

function par_convert(px, conversion_type)

    if conversion_type == 1
        # converts structure array to vector
        p = [px.smpd[1], px.smpd[2], sqrt(px.sdeps)];
            
    elseif conversion_type == 2
        # converts vector to structure array
        p = EstimatedParameters(px[1] * px[1], [px[2], px[3] * px[3]])        
        return p
    else
        println("unvalid conversion type.")
    end

end
