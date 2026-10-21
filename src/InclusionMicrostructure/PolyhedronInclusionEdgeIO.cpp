/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <PolyhedronInclusionEdgeIO.h>

namespace model
{

    PolyhedronInclusionEdgeIO::PolyhedronInclusionEdgeIO(const size_t &inclusionID_in, const size_t &faceID_in,
                                                         const size_t &sourceID_in, const size_t &sinkID_in)
        : inclusionID(inclusionID_in), faceID(faceID_in), sourceID(sourceID_in), sinkID(sinkID_in)
    {

        //            assert(0 && "FINISH HERE, THIS MUST BE COMPATIBLE WITH ID READER");
    }

    PolyhedronInclusionEdgeIO::PolyhedronInclusionEdgeIO() : inclusionID(0), faceID(0), sourceID(0), sinkID(0)
    {
    }

    PolyhedronInclusionEdgeIO::PolyhedronInclusionEdgeIO(std::stringstream &ss)
        : inclusionID(0), faceID(0), sourceID(0), sinkID(0)
    {
        ss >> inclusionID;
        ss >> faceID;
        ss >> sourceID;
        ss >> sinkID;
    }

    void PolyhedronInclusionEdgeIO::write(std::ostream &os) const
    {
        os << inclusionID
           << "\t"
           /**/
           << faceID
           << "\t"
           /**/
           << sourceID
           << "\t"
           /**/
           << sinkID;
        return;
    }
} // namespace model
