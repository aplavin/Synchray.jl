module AdaptExt

using Synchray: CombinedMedium
import Adapt

Adapt.adapt_structure(to, cm::CombinedMedium) = CombinedMedium(Adapt.adapt(to, cm.objects))

end
