import HeroSection from "@/components/HeroSection";
import MetricsSection from "@/components/MetricsSection";
import GeneExpressionChart from "@/components/GeneExpressionChart";
import CellClusteringSection from "@/components/CellClusteringSection";
import SpatialMapSection from "@/components/SpatialMapSection";
import MethodologySection from "@/components/MethodologySection";
import ConclusionSection from "@/components/ConclusionSection";
import Footer from "@/components/Footer";

const Index = () => {
  return (
    <main className="min-h-screen bg-background">
      <HeroSection />
      <MetricsSection />
      <GeneExpressionChart />
      <CellClusteringSection />
      <SpatialMapSection />
      <MethodologySection />
      <ConclusionSection />
      <Footer />
    </main>
  );
};

export default Index;
