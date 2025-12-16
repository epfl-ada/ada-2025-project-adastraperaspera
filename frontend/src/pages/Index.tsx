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
import RQ1Section from "@/components/RQ1Section";
import RQ2Section from "@/components/RQ2Section";
import RQ3Section from "@/components/RQ3Section";

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
      <RQ1Section />
      <RQ2Section />
      <RQ3Section />
      <SpatialMapSection />
      <MethodologySection />
      <ConclusionSection />
      <Footer />
    </main>
  );
};

export default Index;
