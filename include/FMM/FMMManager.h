/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_FMMManager_H_
#define model_FMMManager_H_
#include <map>
#include <memory>
#include <string>
#include <typeindex>
#include <stdexcept>
#include <utility>
namespace model
{
    // Host-only registry. A family contains matching full/short/long kernels.
    // Registration and evaluation are serialized by the owning physics scheduler.
    // Destruction must precede Kokkos finalization if device interactions exist.
    class FMMManager
    {
        std::shared_ptr<void> runtimeOwner_;
        struct Entry
        {
            std::type_index type;
            std::shared_ptr<const void> object;
        };
        std::map<std::string,Entry> kernels_,interactions_;
    public:
        template<class Family> class KernelHandle
        {
            friend class FMMManager;
            std::shared_ptr<const Family> family_;
            explicit KernelHandle(std::shared_ptr<const Family> f):family_(std::move(f)) {}
        public:
            const Family& family() const { return *family_; }
        };
        FMMManager()=default;
        void keepRuntimeAlive(std::shared_ptr<void> owner) { runtimeOwner_=std::move(owner); }
        FMMManager(const FMMManager&)=delete;
        FMMManager& operator=(const FMMManager&)=delete;
        template<class Family> KernelHandle<Family> registerKernel(const std::string& name,Family family)
        {
            if(name.empty() || kernels_.count(name)) { throw std::invalid_argument("Empty or duplicate FMM kernel name: "+name); }
            auto p=std::make_shared<const Family>(std::move(family));
            kernels_.emplace(name,Entry{typeid(Family),p});
            return KernelHandle<Family>(p);
        }
        template<class Family> KernelHandle<Family> kernel(const std::string& name) const
        {
            const auto& e=kernels_.at(name);
            if(e.type!=std::type_index(typeid(Family))) { throw std::invalid_argument("FMM kernel type mismatch: "+name); }
            return KernelHandle<Family>(std::static_pointer_cast<const Family>(e.object));
        }
        // Interaction determines the backend; no Kokkos dependency in this registry.
        template<class Interaction,class Family,class... Args>
        std::shared_ptr<Interaction> registerInteraction(const std::string& name,const KernelHandle<Family>& kernel,Args&&... args)
        {
            if(name.empty() || interactions_.count(name)) { throw std::invalid_argument("Empty or duplicate FMM interaction name: "+name); }
            bool owned=false;
            for(const auto& e:kernels_) { owned=owned || e.second.object.get()==kernel.family_.get(); }
            if(!owned) { throw std::invalid_argument("FMM kernel belongs to another registry"); }
            auto p=std::make_shared<Interaction>(kernel.family_,std::forward<Args>(args)...);
            interactions_.emplace(name,Entry{typeid(Interaction),p});
            return p;
        }
        // Existing caller handles remain valid snapshots; release them to free memory.
        void unregisterInteraction(const std::string& name) { interactions_.erase(name); }
        void unregisterKernel(const std::string& name) { kernels_.erase(name); }
        void clearInteractions() { interactions_.clear(); }
        std::size_t kernelCount() const { return kernels_.size(); }
        std::size_t interactionCount() const { return interactions_.size(); }
    };
}
#endif
