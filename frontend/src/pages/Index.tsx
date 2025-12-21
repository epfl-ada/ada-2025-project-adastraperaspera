import HeroSection from "@/components/HeroSection";
import IntroductionSection from "@/components/IntroductionSection";
import IntroductionSection2 from "@/components/IntroductionSection2";
import MetricsSection from "@/components/MetricsSection";
import GeneExpressionChart from "@/components/GeneExpressionChart";
import GeneExpressionSection2 from "@/components/GeneExpressionChart2";
import CellClusteringSection from "@/components/CellClusteringSection";
import SpatialMapSection from "@/components/SpatialMapSection";
import MethodologySection from "@/components/MethodologySection";
import ConclusionSection from "@/components/ConclusionSection";
import MiscroscopySection from "@/components/MicroscopySection";
import MicroscopySection2 from "@/components/MicroscopySection2";
import PreprocessingSection from "@/components/PreprocessingSection";
import Footer from "@/components/Footer";
import RQSection from "@/components/RQSection";  
import RQ1Section from "@/components/RQ1Section";
import RQ2Section from "@/components/RQ2Section";
import RQ3Section from "@/components/RQ3Section";
import RQ4Section from "@/components/RQ4Section";
import RQ5Section from "@/components/RQ5Section";
import DiscussionSection from "@/components/DiscussionSection";
import DiscussionSection2 from "@/components/DiscussionSection2";
import Chapter1 from "@/components/Chapter1";
import Chapter2 from "@/components/Chapter2";
import Chapter3 from "@/components/Chapter3";
import Chapter4 from "@/components/Chapter4";
import Chapter5 from "@/components/Chapter5";

const Index = () => {
  return (
    <main className="min-h-screen bg-background">
      <HeroSection />
      <IntroductionSection />
      <PreprocessingSection />
      {/*<MetricsSection /> 
      <MicroscopySection/>*/}
      {/*<GeneExpressionChart />
      <GeneExpressionSection2 />
      <CellClusteringSection />*/}
      {/*<RQSection />*/}
      <RQ1Section />
      <RQ2Section />
      <RQ3Section />
      <RQ4Section />
      <RQ5Section />
      {/*<SpatialMapSection />
      <MethodologySection />
      <ConclusionSection />*/}
      {/*<DiscussionSection */}
      <DiscussionSection2/>
      <IntroductionSection2/>
      <Chapter1 />
      <Chapter2 />
      <Chapter3 />
      <Chapter4 />
      <Chapter5 />
      <Footer />
    </main>
  );
};

export default Index;
