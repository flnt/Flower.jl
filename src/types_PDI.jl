#PDI
struct PC_tree_t
    # The tree status
    status::Cint
    # The document containing the tree
    document::Ptr{Cvoid} #Ptr{Cstring} ?
    # the node inside the tree
    node ::Ptr{Cvoid}
end

# These are the values for enum PDI_inout_t 
PDI_NONE  =  Int32(0)
PDI_IN    =  Int32(1)
PDI_OUT   =  Int32(2)
PDI_INOUT =  Int32(3)