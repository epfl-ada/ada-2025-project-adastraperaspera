import HeroSection from "@/components/HeroSection";
import IntroductionSection from "@/components/IntroductionSection";
import MetricsSection from "@/components/MetricsSection";
import GeneExpressionChart from "@/components/GeneExpressionChart";
import CellClusteringSection from "@/components/CellClusteringSection";
import SpatialMapSection from "@/components/SpatialMapSection";
import MethodologySection from "@/components/MethodologySection";
import ConclusionSection from "@/components/ConclusionSection";
import MiscroscopySection from "@/components/MicroscopySection";
import Footer from "@/components/Footer";
import RQSection from "@/components/RQSection";  

const Index = () => {
  return (
    <main className="min-h-screen bg-background">
      <HeroSection />
      <IntroductionSection />
      {/*<MetricsSection /> */}
      <MiscroscopySection/>
      <GeneExpressionChart />
      <CellClusteringSection />
      <RQSection />
      <SpatialMapSection />
      <MethodologySection />
      <ConclusionSection />
      <Footer />
    </main>
  );
};

export default Index;
